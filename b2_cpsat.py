"""CP-SAT model for the interurban Bus Network Design Problem (Roca-Riu et al., 2012).

CP-SAT is an exact solver (clause learning + LP + a portfolio of Large Neighbourhood
Search workers). Run with a time limit it behaves as an anytime heuristic: LNS
workers improve the incumbent, and the best bound gives a (weak) proven gap.
`--lns-only` disables the proof workers and makes it a pure LNS heuristic.

Why a separate model: CP-SAT only accepts integer variables and integer
coefficients (also in the objective in this OR-Tools version), so all data are
scaled to integers. It also offers global constraints that fit this problem:

  s_p   in S_p            stop index at position p = 1..N
  x_pi  <=> (s_p == i)    one-hot channel, exactly one per position
  (s_{p-1}, s_p, t_p, d_p) in Table_p                             (allowed arcs, time, distance)
  sum_p x_pi <= 1                                                 (no stop visited twice)
  A_p = t_1 + ... + t_p                                           (arrival time, exact)
  C_jp = sum_i cost_ji x_pi                                       (stop -> POI perceived time)
  U_j = min_p (W*A_p + C_jp)                                      (lower level, AoN)
  t_1+..+t_N <= Element(b, Tmax)                                  (Eq. 3, buses)
  min rho_O (ct t0 b + cd sum d_p) + rho_U VoT/60 sum_j d_j U_j    (Eqs. 1-7)

"At most N" stops: position p may hold a dummy end value E_p (time/distance 0,
POI cost BIG); once a route ends it stays ended.

Warm start: a greedy route (best first stop, then best extension) is given as a
complete hint, so CP-SAT starts from a feasible solution.

Scaling: times in 0.01 min, distances in metres, objective in 1e-6 EUR. Reported
costs are recomputed exactly with `evaluate_route`; the rounding error is < 0.01%.

Performance notes: the arc tables have 1.6e5 (N=2) to 1.6e6 (N=8) tuples on the
complete 922-stop instance. CP-SAT presolve then takes minutes, so presolve is off
by default. On ~200-450 stops CP-SAT finds the optimum (or near) in 15-30 s; on the
complete instance it stalls far from the optimum: reduce the stop set first.

Usage
-----
python b2_cpsat.py --n-stops 3 5 8 --time-limit 120
python b2_cpsat.py --n-stops 5 --time-limit 60 --lns-only --no-hint
python b2_cpsat.py --n-stops 3 --sample-stops 200 --seed 0 --time-limit 60
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import pandas as pd
from ortools.sat.python import cp_model

from b1_exact_ortools import position_sets
from bndp_common import GATE_ID, ROOT, Instance, Params, evaluate_route, load_instance

TIME_SCALE = 100       # 1 unit = 0.01 min
W_SCALE = 100          # resolution of the w_iu weight
OBJ_SCALE = 1_000_000  # 1 unit = 1e-6 EUR


def greedy_route(inst: Instance, prm: Params, S: list[set[str]], n_stops: int, stops_mode: str) -> list[str] | None:
    """Best first stop, then repeatedly append the stop that most reduces the objective."""
    out = inst.out_arcs()
    route: list[str] = []
    best_obj = math.inf
    for p in range(1, n_stops + 1):
        tail = route[-1] if route else GATE_ID
        cands = [l for l, _, _ in out.get(tail, []) if l in S[p] and l not in route]
        if not cands:
            break
        obj, nxt = min((evaluate_route(inst, prm, route + [l])["objective"], l) for l in cands)
        if stops_mode == "max" and obj >= best_obj:
            break
        route.append(nxt)
        best_obj = obj
    if stops_mode == "eq" and len(route) != n_stops:
        return None
    return route or None


class Progress(cp_model.CpSolverSolutionCallback):
    def __init__(self, s_vars, idx_to_stop, t0, verbose):
        super().__init__()
        self.s_vars, self.idx_to_stop, self.t0, self.verbose = s_vars, idx_to_stop, t0, verbose
        self.trajectory = []

    def on_solution_callback(self):
        obj = self.objective_value / OBJ_SCALE
        route = [self.idx_to_stop[self.value(v)] for v in self.s_vars if self.value(v) in self.idx_to_stop]
        t = time.perf_counter() - self.t0
        self.trajectory.append({"time_s": t, "objective_scaled": obj, "route": route})
        if self.verbose:
            print(f"    {t:8.2f}s  obj={obj:.2f}  bound={self.best_objective_bound / OBJ_SCALE:.2f}  {route}", flush=True)


def build_and_solve(
    inst: Instance,
    prm: Params,
    n_stops: int,
    stops_mode: str = "eq",
    time_limit_s: float = 60.0,
    workers: int = 8,
    seed: int = 0,
    lns_only: bool = False,
    use_hint: bool = True,
    rel_gap: float | None = None,
    log: bool = False,
    verbose: bool = True,
    presolve: str = "off",
) -> dict:
    t_start = time.perf_counter()
    res = {"solver": "CP-SAT" + ("-LNS" if lns_only else ""), "n_stops": n_stops, "stops_mode": stops_mode}

    S = position_sets(inst, n_stops, stops_mode)
    P = list(range(1, n_stops + 1))
    if not S[1] or (stops_mode == "eq" and not S[n_stops]):
        return {**res, "status": "INFEASIBLE", "build_time_s": 0.0, "solve_time_s": 0.0}

    used = sorted(set().union(*S[1:]))
    stop_idx = {k: i for i, k in enumerate(used)}
    idx_to_stop = {i: k for k, i in stop_idx.items()}
    n = len(used)
    end = {p: n + p - 1 for p in P}  # dummy end value of position p ("at most N" mode)
    out = inst.out_arcs()
    ti = lambda t: int(round(t * TIME_SCALE))  # noqa: E731
    di = lambda d: int(round(d * 1000))  # noqa: E731

    m = cp_model.CpModel()

    # --- Positions: s_p (stop index) with a one-hot channel x[p][i] <=> (s_p == i)
    s, x, t, d, tables, t_ub = {}, {}, {}, {}, {}, {}
    for p in P:
        tails = [GATE_ID] if p == 1 else sorted(S[p - 1])
        tup = []
        for k in tails:
            for l, tt, dd in out.get(k, []):
                if l in S[p]:
                    tup.append(((stop_idx[k],) if p > 1 else ()) + (stop_idx[l], ti(tt), di(dd)))
        if stops_mode == "max" and p > 1:
            tup += [(stop_idx[k], end[p], 0, 0) for k in tails]
            if p > 2:
                tup.append((end[p - 1], end[p], 0, 0))
        tables[p] = tup
        vals = sorted({row[-3] for row in tup})
        s[p] = m.new_int_var_from_domain(cp_model.Domain.from_values(vals), f"s{p}")
        x[p] = {}
        for i in vals:
            x[p][i] = m.new_bool_var(f"x{p}_{i}")
            m.add(s[p] == i).only_enforce_if(x[p][i])
            m.add(s[p] != i).only_enforce_if(~x[p][i])
        m.add_exactly_one(x[p].values())
        # --- Arcs between consecutive positions (Eqs. 4, 5, 8-10)
        t_ub[p] = max(r[-2] for r in tup)
        t[p] = m.new_int_var(min(r[-2] for r in tup), t_ub[p], f"t{p}")
        d[p] = m.new_int_var(min(r[-1] for r in tup), max(r[-1] for r in tup), f"d{p}")
        cols = ([s[p - 1]] if p > 1 else []) + [s[p], t[p], d[p]]
        m.add_allowed_assignments(cols, tup)
    for i in range(n):  # each stop visited at most once
        lits = [x[p][i] for p in P if i in x[p]]
        if len(lits) > 1:
            m.add_at_most_one(lits)

    A, A_ub = {}, {}
    for p in P:
        A_ub[p] = sum(t_ub[q] for q in P if q <= p)
        A[p] = m.new_int_var(0, A_ub[p], f"A{p}")
        m.add(A[p] == sum(t[q] for q in P if q <= p))
    t_city, t_city_ub = A[P[-1]], A_ub[P[-1]]

    # --- Buses (Eq. 3): T_C <= (b/f - T_A - T_rest) * 60 min
    b_min = math.ceil(prm.frequency * (prm.t_outside_h + prm.t_rest_h) - 1e-9)
    b_max = math.ceil(prm.frequency * (prm.t_outside_h + prm.t_rest_h + t_city_ub / TIME_SCALE / 60))
    thr = [math.floor((bb / prm.frequency - prm.t_outside_h - prm.t_rest_h) * 60 * TIME_SCALE + 1e-6)
           for bb in range(b_max + 1)]
    b = m.new_int_var(b_min, b_max, "buses")
    b_thr = m.new_int_var(min(thr), max(thr), "bus_time_cap")
    m.add_element(b, thr, b_thr)
    m.add(t_city <= b_thr)

    # --- Lower level (Eq. 6): each POI uses the position minimising W*A_p + C_jp,
    #     with C_jp = sum_i cost_ji * x[p][i] (linear, no element encoding)
    W = int(round(prm.w_iu * W_SCALE))
    real_max = max(c for k in used for c in inst.poi_cost[k].values())
    BIG = int(round(10 * real_max * TIME_SCALE * W_SCALE))
    cost = {}
    U = {}
    hint_terms = []
    for j in inst.pois:
        cost[j] = [BIG] * (n + n_stops)
        for k, i in stop_idx.items():
            c = inst.poi_cost[k].get(j)
            if c is not None:
                cost[j][i] = int(round(c * TIME_SCALE * W_SCALE))
        terms = []
        for p in P:
            c_jp = m.new_int_var(min(cost[j][i] for i in x[p]), max(cost[j][i] for i in x[p]), f"c_{j}_{p}")
            m.add(c_jp == sum(cost[j][i] * lit for i, lit in x[p].items()))
            u_jp = m.new_int_var(0, W * A_ub[p] + BIG, f"u_{j}_{p}")
            m.add(u_jp == W * A[p] + c_jp)
            terms.append(u_jp)
            hint_terms.append((j, p, c_jp, u_jp))
        U[j] = m.new_int_var(0, W * t_city_ub + BIG, f"U_{j}")
        m.add_min_equality(U[j], terms)

    # --- Objective (Eq. 7), scaled to integers
    coef_b = int(round(prm.rho_operator * prm.ct * prm.t0 * OBJ_SCALE))
    coef_d = int(round(prm.rho_operator * prm.cd / 1000 * OBJ_SCALE))
    coef_u = {j: int(round(prm.rho_user * prm.vot / 60 * inst.demand[j] / (TIME_SCALE * W_SCALE) * OBJ_SCALE))
              for j in inst.pois}
    m.minimize(coef_b * b + coef_d * sum(d.values()) + sum(coef_u[j] * U[j] for j in inst.pois))

    # --- Warm start: complete hint (all model variables) from the greedy route
    hint = greedy_route(inst, prm, S, n_stops, stops_mode) if use_hint else None
    if hint:
        sv = {p: stop_idx[hint[p - 1]] if p <= len(hint) else end[p] for p in P}
        prev, acc, arr = GATE_ID, 0, {}
        for p in P:
            m.add_hint(s[p], sv[p])
            for i, lit in x[p].items():
                m.add_hint(lit, i == sv[p])
            if p <= len(hint):
                tt, dd = inst.arcs[(prev, hint[p - 1])]
                prev = hint[p - 1]
            else:
                tt, dd = 0.0, 0.0
            acc += ti(tt)
            m.add_hint(t[p], ti(tt))
            m.add_hint(d[p], di(dd))
            m.add_hint(A[p], acc)
            arr[p] = acc
        bb = next(v for v in range(b_min, b_max + 1) if acc <= thr[v])
        m.add_hint(b, bb)
        m.add_hint(b_thr, thr[bb])
        for j, p, c_jp, u_jp in hint_terms:
            m.add_hint(c_jp, cost[j][sv[p]])
            m.add_hint(u_jp, W * arr[p] + cost[j][sv[p]])
        for j in inst.pois:
            m.add_hint(U[j], min(W * arr[p] + cost[j][sv[p]] for p in P))
        res["hint_route"] = hint
        res["hint_objective"] = evaluate_route(inst, prm, hint)["objective"]

    # --- Solve
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit_s
    solver.parameters.num_workers = workers
    solver.parameters.random_seed = seed
    solver.parameters.log_search_progress = log
    if lns_only:
        solver.parameters.use_lns_only = True
    if rel_gap is not None:
        solver.parameters.relative_gap_limit = rel_gap
    if presolve == "off":
        solver.parameters.cp_model_presolve = False
    elif presolve == "light":
        solver.parameters.max_presolve_iterations = 1
        solver.parameters.cp_model_probing_level = 0
        solver.parameters.symmetry_level = 0

    build_s = time.perf_counter() - t_start
    n_tuples = sum(len(v) for v in tables.values())
    if verbose:
        print(f"[CP-SAT N={n_stops}] {n} stops, {n_tuples} arc tuples, "
              f"{len(m.proto.variables)} vars, {len(m.proto.constraints)} constraints (built in {build_s:.1f}s)"
              + (f", hint obj={res['hint_objective']:.2f}" if hint else ""), flush=True)
    cb = Progress([s[p] for p in P], idx_to_stop, time.perf_counter(), verbose)
    status = solver.solve(m, cb)
    solve_s = time.perf_counter() - t_start - build_s

    res.update({
        "status": solver.status_name(status),
        "build_time_s": build_s,
        "solve_time_s": solve_s,
        "n_vars": len(m.proto.variables),
        "n_constraints": len(m.proto.constraints),
    })
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return res

    route = [idx_to_stop[solver.value(s[p])] for p in P if solver.value(s[p]) in idx_to_stop]
    obj_scaled = solver.objective_value / OBJ_SCALE
    bound = solver.best_objective_bound / OBJ_SCALE
    ev = evaluate_route(inst, prm, route)
    first = cb.trajectory[0]["time_s"] if cb.trajectory else None
    res.update({
        "cpsat_objective": obj_scaled,
        "best_bound": bound,
        "gap": abs(obj_scaled - bound) / max(abs(obj_scaled), 1e-9),
        "time_to_first_s": first,
        "time_to_best_s": cb.trajectory[-1]["time_s"] if cb.trajectory else None,
        "route_names": [inst.names.get(k, k) for k in route],
        **ev,
    })
    if abs(ev["objective"] - obj_scaled) > 1e-3 * abs(ev["objective"]):
        print(f"  WARNING: CP-SAT objective {obj_scaled:.4f} != evaluated objective {ev['objective']:.4f}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stops", type=int, nargs="+", default=[2])
    ap.add_argument("--stops-mode", choices=["eq", "max"], default="eq")
    ap.add_argument("--time-limit", type=float, default=1200, help="seconds per solve")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cp-seed", type=int, default=0, help="CP-SAT random seed")
    ap.add_argument("--lns-only", action="store_true", help="pure LNS heuristic (no proof workers)")
    ap.add_argument("--no-hint", action="store_true", help="do not warm start with the greedy route")
    ap.add_argument("--gap", type=float, default=None, help="stop at this relative gap")
    ap.add_argument("--log", action="store_true", help="CP-SAT search log")
    ap.add_argument("--presolve", choices=["full", "light", "off"], default="off",
                    help="CP-SAT presolve effort; full/light can take minutes on large stop sets")
    ap.add_argument("--out-dir", default=str(ROOT / "results/cpsat"))
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
    name = "CPSAT-LNS" if args.lns_only else "CPSAT"
    summary = []
    for n in args.n_stops:
        res = build_and_solve(inst, prm, n, args.stops_mode, args.time_limit, args.workers, args.cp_seed,
                              args.lns_only, not args.no_hint, args.gap, args.log, presolve=args.presolve)
        res.update({"sample_stops": args.sample_stops, "seed": args.seed, "cp_seed": args.cp_seed,
                    "time_limit_s": args.time_limit})
        (out_dir / f"{name}_N{n}_{args.stops_mode}{tag}_cp{args.cp_seed}.json").write_text(
            json.dumps({"params": vars(prm), **res}, indent=2, ensure_ascii=False))
        summary.append(res)
        if "objective" in res:
            print(f"  {res['status']} obj={res['objective']:.2f} (bound {res['best_bound']:.2f}, gap "
                  f"{100 * res['gap']:.2f}%) C_O={res['operator_cost']:.2f} C_U={res['user_cost']:.2f} "
                  f"buses={res['buses']} T_C={res['city_time_min']:.1f}min best found at "
                  f"{res['time_to_best_s']:.1f}s / {res['solve_time_s']:.1f}s")
            print("  route: Gate -> " + " -> ".join(res["route_names"]))
        else:
            print(f"  {res['status']} after {res['solve_time_s']:.1f}s")

    cols = ["solver", "sample_stops", "seed", "cp_seed", "n_stops", "status", "objective", "best_bound", "gap",
            "operator_cost", "user_cost", "buses", "city_time_min", "length_km", "avg_user_time_min",
            "time_to_first_s", "time_to_best_s", "solve_time_s", "time_limit_s", "route"]
    df = pd.DataFrame([{c: r.get(c) for c in cols} for r in summary])
    df["route"] = df["route"].apply(lambda r: " -> ".join(r) if isinstance(r, list) else r)
    df.to_csv(out_dir / "summary.csv", index=False, mode="a", header=not (out_dir / "summary.csv").exists())
    print(df.drop(columns="route").to_string(index=False))


if __name__ == "__main__":
    main()
