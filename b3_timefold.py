"""Timefold Solver (metaheuristics) for the interurban Bus Network Design Problem (Roca-Riu et al., 2012).

Requirements
------------
Timefold for Python runs on the JVM and only supports Python 3.10-3.12:
  - JDK 17+ (brew install openjdk@21; JAVA_HOME is set automatically if missing)
  - the dedicated venv:  .venv-timefold/bin/python b3_timefold.py ...
    (created with: ~/.pyenv/versions/3.12.3/bin/python -m venv .venv-timefold &&
     .venv-timefold/bin/pip install timefold==1.24.0b0 pandas)

Model
-----
Planning entities: N route positions ("slots"); slot p holds one stop (or END in
at-most-N mode). The route is the sequence of non-END slots in position order,
starting at the gate.

Hard score (feasibility, must reach 0):
  - missing arcs (gate -> first stop, consecutive stops)
  - stops visited twice
  - exactly-N mode: missing stops; at-most-N mode: an empty route
Soft score: -(objective of Eq. 7) in euro cents, with the paper's cost functions
(Eqs. 1-6) and the All-or-Nothing lower level (each POI alights at its best stop).
Reported results are recomputed with `evaluate_route`, as for HiGHS and CP-SAT.

Timefold's generic moves on this model match the paper's Tabu Search moves:
  ChangeMove (slot p -> another stop)  = Mov1 (change terminal) / Mov4 (change intermediate stop)
  ChangeMove (slot p -> END)           = Mov2 (shorten), at-most-N mode
  ChangeMove (END slot -> stop)        = Mov3 (extend), at-most-N mode
  SwapMove (two slots)                 = reorder two stops

Algorithms (construction heuristic FIRST_FIT = greedy stop by stop, then local search):
  SA  Simulated Annealing: accepts worse moves with probability exp(-delta/T); starting
      temperature --sa-temp in EUR (default 200; 50 did better at N=8, 200 at N=3-5). Best performer.
  TS  Tabu Search (the paper's method). Each step evaluates --ts-accepted moves and applies the
      best non-tabu one, even if worse. Tabus: --tabu-entity (a just-changed position is frozen;
      default N-1, needed to leave the greedy local optimum) and --tabu-value (a just-inserted
      stop cannot be inserted again, ~ the paper's mortality H > 18). Timefold has no exact
      equivalent of the paper's "removed stop cannot return" tabu (undo-move tabu is rejected).
  LA, GD, HC  Late Acceptance, Great Deluge, Hill Climbing. They only accept moves at least as
      good as a reference score; the greedy start is a strict local optimum (no equal-cost
      neighbours, ~72% of stop pairs have no arc), so they usually never move. Baselines only.

Times (time_to_best_s, solve_time_s) are on Timefold's clock; overhead_s is the JVM start-up
and Python -> Java conversion (~10-25 s on the complete instance), not part of the search.
Timefold runs single-threaded (multithreaded solving is Enterprise-only).

Usage
-----
.venv-timefold/bin/python b3_timefold.py --n-stops 3 5 8 --time-limit 120           # TS and SA
.venv-timefold/bin/python b3_timefold.py --n-stops 5 --algorithms SA --sa-temp 100 --rs 0 1 2   # 3 seeds
.venv-timefold/bin/python b3_timefold.py --n-stops 3 --sample-stops 200 --seed 0 --time-limit 30
"""

from __future__ import annotations

import os
import subprocess
import sys

if sys.version_info >= (3, 13):
    sys.exit("Timefold supports Python 3.10-3.12. Run: .venv-timefold/bin/python b3_timefold.py ...")


def _ensure_java_home() -> None:
    if os.environ.get("JAVA_HOME"):
        return
    try:
        home = subprocess.run(["/usr/libexec/java_home", "-v", "17+"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        home = "/opt/homebrew/opt/openjdk@21/libexec/openjdk.jdk/Contents/Home"
    if not os.path.isdir(home):
        sys.exit("No JDK 17+ found. Install one (brew install openjdk@21) or set JAVA_HOME.")
    os.environ["JAVA_HOME"] = home


_ensure_java_home()

import argparse  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import time  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Annotated  # noqa: E402

import pandas as pd  # noqa: E402
from timefold.solver import SolverFactory  # noqa: E402
from timefold.solver.config import (  # noqa: E402
    Duration,
    EnvironmentMode,
    ScoreDirectorFactoryConfig,
    SolverConfig,
    TerminationConfig,
)
from timefold.solver.domain import (  # noqa: E402
    PlanningEntityCollectionProperty,
    PlanningId,
    PlanningScore,
    PlanningVariable,
    ProblemFactCollectionProperty,
    ProblemFactProperty,
    ValueRangeProvider,
    planning_entity,
    planning_solution,
)
from timefold.solver.score import (  # noqa: E402
    Constraint,
    ConstraintCollectors,
    ConstraintFactory,
    HardSoftScore,
    constraint_provider,
)

from bndp_common import GATE_ID, ROOT, Instance, Params, evaluate_route, load_instance  # noqa: E402

MISSING_ARC_T = 30.0    # min, only used to keep the soft score defined on infeasible routes
MISSING_ARC_D = 10.0    # km
MISSING_POI = 1.0e6     # min


# --------------------------------------------------------------------------- domain
# Functions run inside the JVM (translated by Timefold), so they cannot read module
# globals: all data lives on the problem facts.

@dataclass
class Stop:
    id: Annotated[int, PlanningId]
    key: str
    gate_t: float               # gate -> stop (min), -1 if no arc
    gate_d: float               # km
    nxt_t: dict                 # stop id -> arc time (min)
    nxt_d: dict                 # stop id -> arc length (km)
    poi: list                   # perceived time stop -> POI j (min)
    is_end: bool = False        # "no stop" marker (at-most-N mode)


@dataclass
class Ctx:
    id: Annotated[int, PlanningId]
    n_stops: int
    exact: bool
    demand: list
    w_iu: float
    rho_o: float
    rho_u: float
    vot: float
    freq: float
    t_out_h: float
    t_rest_h: float
    cd: float
    ct: float
    t0: float

    def route(self, slots) -> list:
        return [s.stop for s in sorted(slots, key=lambda s: s.id) if s.stop is not None and not s.stop.is_end]

    def hard(self, slots) -> int:
        vs = self.route(slots)
        if len(vs) == 0:
            return self.n_stops
        bad = len(vs) - len(set(v.id for v in vs))
        bad += self.n_stops - len(vs) if self.exact else 0
        if vs[0].gate_t < 0:
            bad += 1
        for i in range(1, len(vs)):
            if vs[i].id not in vs[i - 1].nxt_t:
                bad += 1
        return bad

    def soft(self, slots) -> int:
        vs = self.route(slots)
        if len(vs) == 0:
            return 0
        arr = []
        t = vs[0].gate_t if vs[0].gate_t >= 0 else MISSING_ARC_T
        length = vs[0].gate_d if vs[0].gate_t >= 0 else MISSING_ARC_D
        arr.append(t)
        for i in range(1, len(vs)):
            a, b = vs[i - 1], vs[i]
            if b.id in a.nxt_t:
                t += a.nxt_t[b.id]
                length += a.nxt_d[b.id]
            else:
                t += MISSING_ARC_T
                length += MISSING_ARC_D
            arr.append(t)
        user_min = 0.0
        n_poi = len(self.demand)
        best = [MISSING_POI * 10.0] * n_poi
        for p in range(len(vs)):
            a = self.w_iu * arr[p]
            row = vs[p].poi
            for j in range(n_poi):
                c = a + row[j]
                if c < best[j]:
                    best[j] = c
        dem = self.demand
        for j in range(n_poi):
            user_min += dem[j] * best[j]
        buses = math.ceil((self.t_out_h + t / 60.0 + self.t_rest_h) * self.freq - 1e-9)
        c_op = buses * self.ct * self.t0 + self.cd * length
        c_user = self.vot * user_min / 60.0
        return int(round((self.rho_o * c_op + self.rho_u * c_user) * 100.0))


@planning_entity
@dataclass
class Slot:
    id: Annotated[int, PlanningId]
    stop: Annotated[Stop | None, PlanningVariable] = field(default=None)


@planning_solution
@dataclass
class RoutePlan:
    ctx: Annotated[Ctx, ProblemFactProperty]
    stops: Annotated[list[Stop], ProblemFactCollectionProperty, ValueRangeProvider]
    slots: Annotated[list[Slot], PlanningEntityCollectionProperty]
    score: Annotated[HardSoftScore | None, PlanningScore] = field(default=None)


@constraint_provider
def route_constraints(cf: ConstraintFactory) -> list[Constraint]:
    return [
        cf.for_each_including_unassigned(Slot).group_by(ConstraintCollectors.to_list()).join(Ctx)
        .penalize(HardSoftScore.ONE_HARD, lambda slots, ctx: ctx.hard(slots))
        .as_constraint("Infeasible route"),
        cf.for_each_including_unassigned(Slot).group_by(ConstraintCollectors.to_list()).join(Ctx)
        .penalize(HardSoftScore.ONE_SOFT, lambda slots, ctx: ctx.soft(slots))
        .as_constraint("Total cost (cents)"),
    ]


# --------------------------------------------------------------------------- problem / config

def build_problem(inst: Instance, prm: Params, n_stops: int, stops_mode: str) -> RoutePlan:
    idx = {k: i for i, k in enumerate(inst.stops)}
    stops = []
    out = inst.out_arcs()
    for k, i in idx.items():
        g = inst.arcs.get((GATE_ID, k))
        nxt = [(idx[l], t, d) for l, t, d in out.get(k, []) if l in idx]
        stops.append(Stop(
            id=i, key=k,
            gate_t=g[0] if g else -1.0, gate_d=g[1] if g else 0.0,
            nxt_t={l: t for l, t, _ in nxt}, nxt_d={l: d for l, _, d in nxt},
            poi=[inst.poi_cost[k].get(j, MISSING_POI) for j in inst.pois],
        ))
    if stops_mode == "max":
        # an empty position holds END rather than null: Timefold's value tabu fails on null values
        stops.append(Stop(id=len(stops), key="END", gate_t=-1.0, gate_d=0.0, nxt_t={}, nxt_d={},
                          poi=[MISSING_POI] * len(inst.pois), is_end=True))
    ctx = Ctx(
        id=0, n_stops=n_stops, exact=stops_mode == "eq", demand=[inst.demand[j] for j in inst.pois],
        w_iu=prm.w_iu, rho_o=prm.rho_operator, rho_u=prm.rho_user, vot=prm.vot, freq=prm.frequency,
        t_out_h=prm.t_outside_h, t_rest_h=prm.t_rest_h, cd=prm.cd, ct=prm.ct, t0=prm.t0,
    )
    return RoutePlan(ctx=ctx, stops=stops, slots=[Slot(id=p) for p in range(n_stops)])


def local_search_xml(algorithm: str, tabu_entity: int, tabu_value: int, tabu_move: int, ts_accepted: int,
                     sa_temp_eur: float, late_size: int) -> str:
    if algorithm == "LA":
        ls = (f"<acceptor><lateAcceptanceSize>{late_size}</lateAcceptanceSize></acceptor>"
              "<forager><acceptedCountLimit>1</acceptedCountLimit></forager>")
    elif algorithm == "TS":
        # schema order: entity, value, move tabu
        tabu = ""
        if tabu_entity > 0:
            tabu += f"<entityTabuSize>{tabu_entity}</entityTabuSize>"
        if tabu_value > 0:
            tabu += f"<valueTabuSize>{tabu_value}</valueTabuSize>"
        if tabu_move > 0:
            tabu += f"<moveTabuSize>{tabu_move}</moveTabuSize>"
        ls = (f"<acceptor>{tabu}</acceptor>"
              f"<forager><acceptedCountLimit>{ts_accepted}</acceptedCountLimit></forager>")
    elif algorithm == "SA":
        temp = int(round(sa_temp_eur * 100))
        ls = (f"<acceptor><simulatedAnnealingStartingTemperature>0hard/{temp}soft"
              "</simulatedAnnealingStartingTemperature></acceptor>"
              "<forager><acceptedCountLimit>1</acceptedCountLimit></forager>")
    elif algorithm == "GD":
        ls = "<localSearchType>GREAT_DELUGE</localSearchType>"
    elif algorithm == "HC":
        ls = "<localSearchType>HILL_CLIMBING</localSearchType>"
    else:
        raise ValueError(f"unknown algorithm {algorithm}")
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<solver xmlns="https://timefold.ai/xsd/solver">
  <constructionHeuristic><constructionHeuristicType>FIRST_FIT</constructionHeuristicType></constructionHeuristic>
  <localSearch>{ls}</localSearch>
</solver>"""


ALGO_NAMES = {"LA": "Late Acceptance", "TS": "Tabu Search", "SA": "Simulated Annealing",
              "GD": "Great Deluge", "HC": "Hill Climbing"}


def solve(
    inst: Instance,
    prm: Params,
    n_stops: int,
    algorithm: str = "TS",
    stops_mode: str = "eq",
    time_limit_s: float = 60.0,
    unimproved_s: float | None = None,
    random_seed: int = 0,
    tabu_entity: int = -1,
    tabu_value: int = 20,
    tabu_move: int = 0,
    ts_accepted: int = 500,
    sa_temp_eur: float = 200.0,
    late_size: int = 400,
    verbose: bool = True,
) -> dict:
    t_start = time.perf_counter()
    problem = build_problem(inst, prm, n_stops, stops_mode)
    termination = TerminationConfig(spent_limit=Duration(milliseconds=int(time_limit_s * 1000)))
    if unimproved_s:
        termination = TerminationConfig(
            spent_limit=Duration(milliseconds=int(time_limit_s * 1000)),
            unimproved_spent_limit=Duration(milliseconds=int(unimproved_s * 1000)),
        )
    config = SolverConfig(
        xml_source_text=local_search_xml(algorithm, n_stops - 1 if tabu_entity < 0 else min(tabu_entity, n_stops - 1), tabu_value, tabu_move, ts_accepted,
                                         sa_temp_eur, late_size),
        solution_class=RoutePlan,
        entity_class_list=[Slot],
        score_director_factory_config=ScoreDirectorFactoryConfig(constraint_provider_function=route_constraints),
        termination_config=termination,
        environment_mode=EnvironmentMode.NO_ASSERT,  # reproducible for a given seed, without debug asserts
        random_seed=random_seed,
    )
    solver = SolverFactory.create(config).build_solver()

    # times are taken from Timefold's clock (from solving start, i.e. excluding the
    # JVM start-up and the Python -> Java conversion of the problem: see overhead_s)
    improvements = []

    def on_best(event) -> None:
        sc = event.new_best_score
        t = event.time_spent.total_seconds()
        improvements.append((t, sc.hard_score, sc.soft_score))
        if verbose and sc.hard_score == 0:
            print(f"    {t:8.2f}s  obj={-sc.soft_score / 100:.2f}", flush=True)

    solver.add_event_listener(on_best)
    build_s = time.perf_counter() - t_start
    if verbose:
        print(f"[Timefold {algorithm} N={n_stops}] {len(inst.stops)} stops (built in {build_s:.1f}s)", flush=True)
    t0 = time.perf_counter()
    best = solver.solve(problem)
    wall_s = time.perf_counter() - t0
    try:
        solve_s = solver._delegate.getTimeMillisSpent() / 1000  # Timefold's own clock
    except Exception:
        solve_s = min(time_limit_s, wall_s)

    route = [s.stop.key for s in sorted(best.slots, key=lambda s: s.id) if s.stop is not None and not s.stop.is_end]
    feasible = best.score.hard_score == 0
    res = {
        "solver": f"Timefold-{algorithm}",
        "algorithm": ALGO_NAMES[algorithm],
        "n_stops": n_stops,
        "stops_mode": stops_mode,
        "status": "FEASIBLE" if feasible else "INFEASIBLE",
        "timefold_score": str(best.score),
        "build_time_s": build_s,
        "solve_time_s": solve_s,
        "overhead_s": wall_s - solve_s,
        "random_seed": random_seed,
    }
    feas = [(t, s) for t, h, s in improvements if h == 0]
    if feasible:
        ev = evaluate_route(inst, prm, route)
        res.update({
            "timefold_objective": -best.score.soft_score / 100,
            "time_to_first_s": feas[0][0] if feas else None,
            "time_to_best_s": feas[-1][0] if feas else None,
            "route_names": [inst.names.get(k, k) for k in route],
            **ev,
        })
        if abs(ev["objective"] - res["timefold_objective"]) > 0.02:
            print(f"  WARNING: Timefold objective {res['timefold_objective']:.2f} != evaluated {ev['objective']:.2f}")
    else:
        res["route"] = route
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stops", type=int, nargs="+", default=[2])
    ap.add_argument("--stops-mode", choices=["eq", "max"], default="eq")
    ap.add_argument("--algorithms", nargs="+", default=["TS", "SA"], choices=list(ALGO_NAMES))
    ap.add_argument("--time-limit", type=float, default=3600, help="seconds per solve")
    ap.add_argument("--unimproved", type=float, default=120, help="stop after this many seconds without improvement")
    ap.add_argument("--rs", type=int, nargs="+", default=[0], help="Timefold random seed(s)")
    ap.add_argument("--tabu-entity", type=int, default=-1,
                    help="TS: steps a just-changed route position cannot change again (-1 = N-1, capped at N-1)")
    ap.add_argument("--tabu-value", type=int, default=20,
                    help="TS: steps a just-inserted stop cannot be inserted again (paper: H > 18)")
    ap.add_argument("--tabu-move", type=int, default=0, help="TS: steps an identical move is forbidden")
    ap.add_argument("--ts-accepted", type=int, default=500,
                    help="TS: moves evaluated per step (best non-tabu one is applied)")
    ap.add_argument("--sa-temp", type=float, default=200.0, help="SA: starting temperature in euros")
    ap.add_argument("--late-size", type=int, default=400, help="LA: late acceptance size")
    ap.add_argument("--out-dir", default=str(ROOT / "results/timefold"))
    ap.add_argument("--sample-stops", type=int, default=None, help="random subset of candidate stops (testing)")
    ap.add_argument("--seed", type=int, default=0, help="seed of the random stop subset")
    prm_default = Params()
    for name, val in vars(prm_default).items():
        ap.add_argument(f"--{name.replace('_', '-')}", type=type(val), default=val)
    args = ap.parse_args()

    prm = Params(**{k: getattr(args, k) for k in vars(prm_default)})
    print("Loading data...", flush=True)
    inst = load_instance(sample_stops=args.sample_stops, seed=args.seed)
    print(f"{len(inst.stops)} candidate stops, {len(inst.arcs)} arcs, {len(inst.pois)} POIs", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"_S{args.sample_stops}_seed{args.seed}" if args.sample_stops else ""
    summary = []
    for n in args.n_stops:
        for algo in args.algorithms:
            for rs in args.rs:
                res = solve(inst, prm, n, algo, args.stops_mode, args.time_limit, args.unimproved, rs,
                            args.tabu_entity, args.tabu_value, args.tabu_move, args.ts_accepted, args.sa_temp,
                            args.late_size)
                res.update({"sample_stops": args.sample_stops, "seed": args.seed, "time_limit_s": args.time_limit})
                (out_dir / f"Timefold-{algo}_N{n}_{args.stops_mode}{tag}_rs{rs}.json").write_text(
                    json.dumps({"params": vars(prm), **res}, indent=2, ensure_ascii=False))
                summary.append(res)
                if "objective" in res:
                    print(f"  {res['status']} obj={res['objective']:.2f} C_O={res['operator_cost']:.2f} "
                          f"C_U={res['user_cost']:.2f} buses={res['buses']} T_C={res['city_time_min']:.1f}min "
                          f"best found at {res['time_to_best_s']:.1f}s / {res['solve_time_s']:.1f}s")
                    print("  route: Gate -> " + " -> ".join(res["route_names"]))
                else:
                    print(f"  {res['status']} ({res['timefold_score']}) after {res['solve_time_s']:.1f}s")

    cols = ["solver", "sample_stops", "seed", "random_seed", "n_stops", "status", "objective", "operator_cost",
            "user_cost", "buses", "city_time_min", "length_km", "avg_user_time_min", "time_to_first_s",
            "time_to_best_s", "solve_time_s", "overhead_s", "time_limit_s", "route"]
    df = pd.DataFrame([{c: r.get(c) for c in cols} for r in summary])
    df["route"] = df["route"].apply(lambda r: " -> ".join(r) if isinstance(r, list) else r)
    df.to_csv(out_dir / "summary.csv", index=False, mode="a", header=not (out_dir / "summary.csv").exists())
    print(df.drop(columns="route").to_string(index=False))


if __name__ == "__main__":
    main()
