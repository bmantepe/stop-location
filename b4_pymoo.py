"""pymoo (evolutionary algorithms) for the interurban Bus Network Design Problem (Roca-Riu et al., 2012).

pymoo is a library of population-based metaheuristics: it evolves a population of
routes with selection, crossover and mutation. It is a heuristic (no optimality proof,
no bound). Two algorithms are provided:

  NSGA2  bi-objective (C_O, C_U): returns the Pareto front of operator vs user cost.
         The reported route is the front member with the best Eq. (7) value (a weighted-sum
         optimum is always on the Pareto front), so it is comparable with the other solvers;
         the whole front is stored in the JSON ("pareto_front"). Default and best performer:
         keeping cheap (few buses) and user-friendly routes in the population preserves
         diversity across bus counts, where the GA tends to converge to one bus count.
  GA     single-objective genetic algorithm on Eq. (7): rho_O*C_O + rho_U*C_U

Encoding: a route is an integer vector of length N with the stop index at each position
(-1 = no stop, at-most-N mode only). Random vectors are almost never valid routes (~72%
of stop pairs have no arc), so custom operators keep every individual feasible:

  Sampling   random walks from the gate along existing arcs (+ the greedy route with --greedy)
  Crossover  splice: head of parent A + tail of parent B where the joining arc exists, or
             joined through one bridging stop k (A[i-1] -> k -> B[j+1])
  Mutation   the paper's Tabu Search moves:
             change a stop (Mov1 terminal / Mov4 intermediate), delete a stop and append a new
             terminal (Mov2 + Mov3, keeps N), and in at-most-N mode also shorten (Mov2),
             extend (Mov3) and insert.

Every individual is scored with the paper's cost functions (Eqs. 1-7, All-or-Nothing lower
level); reported results are recomputed with `evaluate_route`, as in the other scripts.

Stopping rules (the run ends at the first one reached; stop_reason is saved):
  --time-limit S   wall-clock cap (default 60; 0 = none)
  --n-gen G        max generations      --n-evals E   max route evaluations
  --stall-gens G   G generations without improving the best Eq. (7) value
  --stall-time S   S seconds without improvement
  --target V       objective <= V (e.g. a known optimum, to measure time-to-target)
Time-based stopping depends on machine speed, so the same seed can give slightly different
results; --n-gen/--n-evals (with --time-limit 0) give exactly reproducible runs.

Usage
-----
python b4_pymoo.py --n-stops 3 5 8 --time-limit 120                 # NSGA2
python b4_pymoo.py --n-stops 5 --algorithms NSGA2 GA --rs 0 1 2
python b4_pymoo.py --n-stops 3 --sample-stops 200 --seed 0 --time-limit 30
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.core.callback import Callback
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.population import Population
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.sampling import Sampling
from pymoo.core.termination import TerminateIfAny
from pymoo.optimize import minimize
from pymoo.termination import get_termination

from b1_exact_ortools import position_sets
from b2_cpsat import greedy_route
from bndp_common import GATE_ID, ROOT, Instance, Params, evaluate_route, load_instance

GATE = -1


class Network:
    """Index-based view of the instance with a fast (vectorised) route cost."""

    def __init__(self, inst: Instance, prm: Params, n_stops: int, stops_mode: str):
        self.inst, self.prm, self.n, self.exact = inst, prm, n_stops, stops_mode == "eq"
        self.keys = list(inst.stops)
        idx = {k: i for i, k in enumerate(self.keys)}
        self.arc: dict[tuple[int, int], tuple[float, float]] = {}
        self.out: dict[int, set[int]] = {GATE: set()}
        self.inn: dict[int, set[int]] = {}
        for (k, l), td in inst.arcs.items():
            a = GATE if k == GATE_ID else idx[k]
            b = idx[l]
            self.arc[(a, b)] = td
            self.out.setdefault(a, set()).add(b)
            self.inn.setdefault(b, set()).add(a)
        S = position_sets(inst, n_stops, stops_mode)
        self.S = [{idx[k] for k in s} for s in S]  # S[p]: stops allowed at position p (1-based)
        self.C = np.full((len(self.keys), len(inst.pois)), np.inf)
        for k, i in idx.items():
            for jj, j in enumerate(inst.pois):
                c = inst.poi_cost[k].get(j)
                if c is not None:
                    self.C[i, jj] = c
        self.demand = np.array([inst.demand[j] for j in inst.pois])
        g = greedy_route(inst, prm, S, n_stops, stops_mode)
        self.greedy = [idx[k] for k in g] if g else None

    def succ(self, a: int) -> set[int]:
        return self.out.get(a, set())

    def pred(self, b: int) -> set[int]:
        return self.inn.get(b, set())

    def valid(self, r: list[int]) -> bool:
        if not r or len(r) > self.n or (self.exact and len(r) != self.n) or len(set(r)) != len(r):
            return False
        return all((a, b) in self.arc for a, b in zip([GATE] + r[:-1], r))

    def cost(self, r: list[int]) -> tuple[float, float, float]:
        """(objective, C_O, C_U) of a valid route, same formulas as evaluate_route."""
        prm = self.prm
        t = length = 0.0
        arr = np.empty(len(r))
        for p, (a, b) in enumerate(zip([GATE] + r[:-1], r)):
            tt, dd = self.arc[(a, b)]
            t += tt
            length += dd
            arr[p] = t
        best = (prm.w_iu * arr[:, None] + self.C[r]).min(axis=0)
        c_user = prm.vot * float(self.demand @ best) / 60
        buses = math.ceil((prm.t_outside_h + t / 60 + prm.t_rest_h) * prm.frequency - 1e-9)
        c_op = buses * prm.ct * prm.t0 + prm.cd * length
        return prm.rho_operator * c_op + prm.rho_user * c_user, c_op, c_user

    def random_route(self, rng: np.random.Generator) -> list[int] | None:
        for _ in range(200):
            length = self.n if self.exact else int(rng.integers(1, self.n + 1))
            r, prev = [], GATE
            for p in range(1, length + 1):
                cand = [k for k in self.succ(prev) if k in self.S[p] and k not in r]
                if not cand:
                    break
                prev = cand[int(rng.integers(len(cand)))]
                r.append(prev)
            if len(r) == length and self.valid(r):
                return r
        return None

    def encode(self, r: list[int]) -> np.ndarray:
        x = np.full(self.n, -1, dtype=int)
        x[: len(r)] = r
        return x

    @staticmethod
    def decode(x) -> list[int]:
        return [int(v) for v in x if v >= 0]


# --------------------------------------------------------------------------- pymoo operators

class RouteProblem(ElementwiseProblem):
    def __init__(self, net: Network, n_obj: int):
        super().__init__(n_var=net.n, n_obj=n_obj, xl=-1, xu=len(net.keys) - 1, vtype=int)
        self.net = net

    def _evaluate(self, x, out, *args, **kwargs):
        r = self.net.decode(x)
        if not self.net.valid(r):  # operators keep routes valid; this is a safety net
            out["F"] = [1e12] * self.n_obj
            return
        obj, c_op, c_user = self.net.cost(r)
        out["F"] = [obj] if self.n_obj == 1 else [c_op, c_user]


class RouteSampling(Sampling):
    def __init__(self, use_greedy: bool):
        super().__init__()
        self.use_greedy = use_greedy

    def _do(self, problem, n_samples, *args, random_state=None, **kwargs):
        net = problem.net
        X = []
        if self.use_greedy and net.greedy:
            X.append(net.encode(net.greedy))
        while len(X) < n_samples:
            r = net.random_route(random_state)
            if r is None:
                raise RuntimeError("could not build a feasible random route")
            X.append(net.encode(r))
        return np.array(X)


class SpliceCrossover(Crossover):
    def __init__(self, prob: float = 0.9):
        super().__init__(n_parents=2, n_offsprings=2, prob=prob)

    @staticmethod
    def children(net: Network, a: list[int], b: list[int]) -> list[list[int]]:
        out = []
        for i in range(1, len(a) + 1):          # keep a[:i]
            head, last = a[:i], a[i - 1]
            hs = set(head)
            js = [i + len(b) - net.n] if net.exact else range(len(b) + 1)
            for j in js:                         # continue with b[j:]
                if j < 0 or j > len(b):
                    continue
                tail = b[j:]
                if len(head) + len(tail) > net.n or (not tail and net.exact):
                    continue
                if not tail:
                    if head != a:
                        out.append(head)
                    continue
                if hs.isdisjoint(tail) and (last, tail[0]) in net.arc:
                    out.append(head + tail)
                if len(tail) >= 2 and hs.isdisjoint(tail[1:]):   # bridge: last -> k -> tail[1]
                    ks = (net.succ(last) & net.pred(tail[1])) - hs - set(tail[1:])
                    if ks:
                        out.append(head + [min(ks)] + tail[1:])
        return [c for c in out if c != a and c != b and net.valid(c)]

    def _do(self, problem, X, *args, random_state=None, **kwargs):
        net = problem.net
        _, n_matings, n_var = X.shape
        Y = np.array(X, copy=True)
        for m in range(n_matings):
            a, b = net.decode(X[0, m]), net.decode(X[1, m])
            for o, (p, q) in enumerate(((a, b), (b, a))):
                cands = self.children(net, p, q)
                if cands:
                    Y[o, m] = net.encode(cands[int(random_state.integers(len(cands)))])
                else:  # no valid splice (e.g. identical parents): mutate instead of copying
                    Y[o, m] = net.encode(mutate_route(net, p, random_state))
        return Y


def mutate_route(net: Network, r: list[int], rng) -> list[int]:
    """One random feasible move of the paper's Tabu Search (Mov1-Mov4)."""
    L = len(r)
    moves = ["change", "delete_append"]
    if not net.exact:
        moves += ["append", "drop_last", "insert", "delete"]
    for _ in range(30):
        mv = moves[int(rng.integers(len(moves)))]
        p = int(rng.integers(L))
        prev = r[p - 1] if p > 0 else GATE
        used = set(r)
        if mv == "change":                       # Mov1 (p = last) / Mov4
            ks = net.succ(prev) - used
            if p < L - 1:
                ks &= net.pred(r[p + 1])
            new = None if not ks else r[:p] + [sorted(ks)[int(rng.integers(len(ks)))]] + r[p + 1:]
        elif mv == "delete_append":              # Mov2 + Mov3, keeps the length
            rest = r[:p] + r[p + 1:]
            if not rest or (p < L - 1 and (prev, r[p + 1]) not in net.arc):
                continue
            ks = net.succ(rest[-1]) - set(rest) - {r[p]}
            new = None if not ks else rest + [sorted(ks)[int(rng.integers(len(ks)))]]
        elif mv == "append" and L < net.n:       # Mov3
            ks = net.succ(r[-1]) - used
            new = None if not ks else r + [sorted(ks)[int(rng.integers(len(ks)))]]
        elif mv == "drop_last" and L > 1:        # Mov2
            new = r[:-1]
        elif mv == "insert" and L < net.n:
            ks = (net.succ(prev) & net.pred(r[p])) - used
            new = None if not ks else r[:p] + [sorted(ks)[int(rng.integers(len(ks)))]] + r[p:]
        elif mv == "delete" and L > 1:
            new = r[:p] + r[p + 1:] if p == L - 1 or (prev, r[p + 1]) in net.arc else None
        else:
            continue
        if new and new != r and net.valid(new):
            return new
    return r


class RouteMutation(Mutation):
    """Applies one move of mutate_route to each offspring (with probability p_mut)."""

    def _do(self, problem, X, *args, random_state=None, **kwargs):
        net = problem.net
        Y = np.array(X, copy=True)
        for i in range(len(X)):
            Y[i] = net.encode(mutate_route(net, net.decode(X[i]), random_state))
        return Y


class Progress(Callback):
    """Tracks the best Eq. (7) value and applies the stall / target stopping rules."""

    def __init__(self, net: Network, n_obj: int, verbose: bool, stall_gens: int | None = None,
                 stall_s: float | None = None, target: float | None = None):
        super().__init__()
        self.net, self.n_obj, self.verbose = net, n_obj, verbose
        self.stall_gens, self.stall_s, self.target = stall_gens, stall_s, target
        self.t0 = time.perf_counter()
        self.best = math.inf
        self.time_to_first = self.time_to_best = None
        self.n_gen = 0
        self.gen_best = 0
        self.stop_reason = None

    def notify(self, algorithm):
        self.n_gen = algorithm.n_gen
        if self.n_obj == 1:
            obj = float(algorithm.opt[0].F[0])
        else:  # best Eq. (7) value in the current non-dominated set
            obj = min(self.net.cost(self.net.decode(x))[0] for x in algorithm.opt.get("X"))
        t = time.perf_counter() - self.t0
        if self.time_to_first is None:
            self.time_to_first = t
        if obj < self.best - 1e-9:
            self.best, self.time_to_best, self.gen_best = obj, t, algorithm.n_gen
            if self.verbose:
                print(f"    {t:8.2f}s  gen {algorithm.n_gen:5d}  obj={obj:.2f}", flush=True)
        if self.target is not None and self.best <= self.target:
            self.stop_reason = "target"
        elif self.stall_gens and algorithm.n_gen - self.gen_best >= self.stall_gens:
            self.stop_reason = "stall_generations"
        elif self.stall_s and t - self.time_to_best >= self.stall_s:
            self.stop_reason = "stall_time"
        if self.stop_reason:
            algorithm.termination.force_termination = True


# --------------------------------------------------------------------------- driver

def _limit_reason(gens, evals, t, n_gen, n_evals, time_limit_s) -> str:
    if n_gen and gens >= n_gen:
        return "n_gen"
    if n_evals and evals >= n_evals:
        return "n_evals"
    return "time_limit" if time_limit_s and t >= time_limit_s - 0.5 else "other"


def add_random_immigrants(algo, net: Network) -> list[int]:
    """If mating cannot create any new unique offspring, pymoo stops the run early; inject
    fresh random routes instead so the run always uses its full budget."""
    orig = algo._infill
    count = [0]

    def infill():
        off = orig()
        if off is None or len(off) == 0:
            algo.termination.force_termination = False
            count[0] += 1
            X = np.array([net.encode(net.random_route(algo.random_state)) for _ in range(algo.n_offsprings)])
            off = Population.new(X=X)
        return off

    algo._infill = infill
    return count


def solve(
    inst: Instance,
    prm: Params,
    n_stops: int,
    algorithm: str = "NSGA2",
    stops_mode: str = "eq",
    time_limit_s: float | None = 60.0,
    random_seed: int = 0,
    n_gen: int | None = None,
    n_evals: int | None = None,
    stall_gens: int | None = None,
    stall_s: float | None = None,
    target: float | None = None,
    pop_size: int = 100,
    p_cross: float = 0.9,
    p_mut: float = 0.5,
    use_greedy: bool = False,
    verbose: bool = True,
) -> dict:
    t_start = time.perf_counter()
    res = {"solver": f"pymoo-{algorithm}", "n_stops": n_stops, "stops_mode": stops_mode,
           "random_seed": random_seed, "pop_size": pop_size, "p_cross": p_cross, "p_mut": p_mut,
           "greedy_start": use_greedy}
    net = Network(inst, prm, n_stops, stops_mode)
    if not net.S[1] or (net.exact and not net.S[n_stops]):
        return {**res, "status": "INFEASIBLE", "build_time_s": 0.0, "solve_time_s": 0.0}

    n_obj = 1 if algorithm == "GA" else 2
    problem = RouteProblem(net, n_obj)
    ops = dict(pop_size=pop_size, sampling=RouteSampling(use_greedy), crossover=SpliceCrossover(p_cross),
               mutation=RouteMutation(prob=p_mut), eliminate_duplicates=True)
    algo = GA(**ops) if algorithm == "GA" else NSGA2(**ops)
    immigrants = add_random_immigrants(algo, net)
    # stop at whichever limit comes first; stall/target rules are applied in Progress
    limits = []
    if time_limit_s:
        limits.append(get_termination("time", time_limit_s))
    if n_gen:
        limits.append(get_termination("n_gen", n_gen))
    if n_evals:
        limits.append(get_termination("n_evals", n_evals))
    if not limits:
        raise ValueError("set at least one of time limit, n_gen or n_evals (stall rules alone may never stop)")
    termination = TerminateIfAny(*limits)
    build_s = time.perf_counter() - t_start
    if verbose:
        print(f"[pymoo {algorithm} N={n_stops}] {len(net.keys)} stops, pop {pop_size}"
              + (f", greedy obj={net.cost(net.greedy)[0]:.2f}" if use_greedy and net.greedy else "")
              + f" (built in {build_s:.1f}s)", flush=True)

    cb = Progress(net, n_obj, verbose, stall_gens, stall_s, target)
    out = minimize(problem, algo, termination, seed=random_seed, callback=cb, verbose=False,
                   copy_algorithm=False)  # keep the immigrant wrapper on the algorithm that runs
    solve_s = time.perf_counter() - cb.t0

    X = np.atleast_2d(out.X)
    routes = [net.decode(x) for x in X]
    best = min(routes, key=lambda r: net.cost(r)[0])
    route = [net.keys[i] for i in best]
    ev = evaluate_route(inst, prm, route)
    res.update({
        "status": "FEASIBLE",
        "build_time_s": build_s,
        "solve_time_s": solve_s,
        "n_gen": cb.n_gen,
        "n_evals": int(algo.evaluator.n_eval),
        "stop_reason": cb.stop_reason or _limit_reason(cb.n_gen, int(algo.evaluator.n_eval), solve_s,
                                                       n_gen, n_evals, time_limit_s),
        "immigrant_injections": immigrants[0],
        "time_to_first_s": cb.time_to_first,
        "time_to_best_s": cb.time_to_best,
        "route_names": [inst.names.get(k, k) for k in route],
        **ev,
    })
    if n_obj == 2:
        front = []
        for r in sorted(routes, key=lambda r: net.cost(r)[1]):
            e = evaluate_route(inst, prm, [net.keys[i] for i in r])
            front.append({k: e[k] for k in ("route", "operator_cost", "user_cost", "objective", "buses",
                                            "city_time_min", "length_km")})
        res["pareto_front"] = front
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stops", type=int, nargs="+", default=[2])
    ap.add_argument("--stops-mode", choices=["eq", "max"], default="eq")
    ap.add_argument("--algorithms", nargs="+", default=["NSGA2"], choices=["GA", "NSGA2"])
    ap.add_argument("--time-limit", type=float, default=60.0,
                    help="max seconds per run (0 = no time cap; then set --n-gen or --n-evals)")
    ap.add_argument("--n-gen", type=int, default=None,
                    help="max generations (exactly reproducible when the time cap is not reached)")
    ap.add_argument("--n-evals", type=int, default=None, help="max route evaluations")
    ap.add_argument("--stall-gens", type=int, default=None, help="stop after this many generations without improvement")
    ap.add_argument("--stall-time", type=float, default=None, help="stop after this many seconds without improvement")
    ap.add_argument("--target", type=float, default=None, help="stop once the objective is <= this value")
    ap.add_argument("--rs", type=int, nargs="+", default=[0], help="random seed(s) of the algorithm")
    ap.add_argument("--pop-size", type=int, default=100)
    ap.add_argument("--p-cross", type=float, default=0.9, help="crossover probability")
    ap.add_argument("--p-mut", type=float, default=0.5, help="mutation probability per offspring")
    ap.add_argument("--greedy", action="store_true",
                    help="put the greedy route in the initial population (tends to take over the GA population)")
    ap.add_argument("--out-dir", default=str(ROOT / "results/pymoo"))
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
                res = solve(inst, prm, n, algo, args.stops_mode, args.time_limit or None, rs, args.n_gen, args.n_evals,
                            args.stall_gens, args.stall_time, args.target, args.pop_size,
                            args.p_cross, args.p_mut, args.greedy)
                res.update({"sample_stops": args.sample_stops, "seed": args.seed,
                            "time_limit_s": args.time_limit or None, "max_gen": args.n_gen,
                            "max_evals": args.n_evals, "stall_gens": args.stall_gens,
                            "stall_time_s": args.stall_time, "target": args.target})
                (out_dir / f"pymoo-{algo}_N{n}_{args.stops_mode}{tag}_rs{rs}.json").write_text(
                    json.dumps({"params": vars(prm), **res}, indent=2, ensure_ascii=False))
                summary.append(res)
                if "objective" in res:
                    print(f"  {res['status']} obj={res['objective']:.2f} C_O={res['operator_cost']:.2f} "
                          f"C_U={res['user_cost']:.2f} buses={res['buses']} T_C={res['city_time_min']:.1f}min "
                          f"best found at {res['time_to_best_s']:.1f}s / {res['solve_time_s']:.1f}s "
                          f"({res['n_gen']} generations, stopped by {res['stop_reason']})"
                          + (f", Pareto front: {len(res['pareto_front'])} routes" if "pareto_front" in res else ""))
                    print("  route: Gate -> " + " -> ".join(res["route_names"]))
                else:
                    print(f"  {res['status']}")

    cols = ["solver", "sample_stops", "seed", "random_seed", "n_stops", "stops_mode", "status", "objective",
            "operator_cost", "user_cost", "buses", "city_time_min", "length_km", "avg_user_time_min",
            "time_to_first_s", "time_to_best_s", "solve_time_s", "stop_reason", "time_limit_s", "max_gen",
            "max_evals", "stall_gens", "stall_time_s", "target", "n_gen", "n_evals", "pop_size",
            "p_cross", "p_mut", "greedy_start", "route"]
    df = pd.DataFrame([{c: r.get(c) for c in cols} for r in summary])
    df["route"] = df["route"].apply(lambda r: " -> ".join(r) if isinstance(r, list) else r)
    df.to_csv(out_dir / "summary.csv", index=False, mode="a", header=not (out_dir / "summary.csv").exists())
    print(df.drop(columns="route").to_string(index=False))


if __name__ == "__main__":
    main()
