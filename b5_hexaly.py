"""Hexaly Optimizer for the interurban Bus Network Design Problem (Roca-Riu et al., 2012).

Hexaly is a commercial solver driven by local search (plus internal bounding). It is used here
as a heuristic: it returns the best route found within the time limit and a lower bound that
is usually weak.

Licence (needed to run anything)
--------------------------------
Download `license.dat` from your account on hexaly.com (Account > Licenses) and either
  - put it in this project folder (it is git-ignored), or
  - pass --license /path/to/license.dat, or
  - set HX_LICENSE_PATH=/path/to/license.dat, or
  - copy it to /opt/hexaly_15_0/license.dat (Hexaly's default location), or
  - set HX_LICENSE_CONTENT to the licence content (the whole text of license.dat).
The Python package itself is `pip install hexaly` (already in the default Python 3.13).

Model
-----
route = list(n)                  ordered list of stop indices (Hexaly list variable)
count(route) == N  (or 1..N)
arrival[i] = gate_t[route[0]] if i == 0 else arrival[i-1] + T[route[i-1]][route[i]]   (recursive array)
T_C = arrival[c-1];  length = gate_d[route[0]] + sum D[route[i-1]][route[i]]           (Eqs. 4-5)
buses = ceil((T_A + T_C/60 + T_rest) * f)                                               (Eq. 3)
best_j = min over positions p of  w_iu*arrival[p] + C[route[p]][j]       (lower level, AoN)
cost   = rho_O*(ct*t0*buses + cd*length) + rho_U*VoT/60*sum_j d_j best_j                (Eqs. 1-7)

Objectives (lexicographic, as in Hexaly's routing templates):
  1) number of missing arcs (gate -> first stop, consecutive stops): 0 = feasible route
  2) cost (Eq. 7)
Treating missing arcs as a first-priority objective rather than a hard constraint lets the
local search move through infeasible routes; only routes with 0 missing arcs are reported.
Reported results are recomputed with `evaluate_route`, as in the other scripts.

Stopping rules (the run ends at the first one reached; stop_reason is saved):
  --time-limit S   wall-clock cap in whole seconds (default 60; 0 = none)
  --iterations I   max Hexaly iterations (moves tried; reproducible with --threads 1)
  --stall-time S   S seconds without improvement (checked every second)
  --target V       objective <= V

Usage
-----
python b5_hexaly.py --n-stops 3 5 8 --time-limit 120
python b5_hexaly.py --n-stops 5 --rs 0 1 2 --threads 8
python b5_hexaly.py --n-stops 3 --sample-stops 60 --seed 2 --time-limit 10 --check   # vs brute force
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import pandas as pd

import hexaly.optimizer as hx
from b1_exact_ortools import position_sets
from b2_cpsat import greedy_route
from bndp_common import GATE_ID, ROOT, Instance, Params, evaluate_route, load_instance

MISSING_ARC_T = 60.0  # min; only steers the search while arcs are missing (objective 1 > 0)
MISSING_ARC_D = 20.0  # km
MISSING_POI = 1.0e5   # min


def set_license(path: str | None) -> str:
    """Point Hexaly to the licence file; returns the path used."""
    if os.environ.get("HX_LICENSE_CONTENT") and not path:
        return "HX_LICENSE_CONTENT environment variable"  # read by Hexaly itself, has priority
    candidates = [path, os.environ.get("HX_LICENSE_PATH"), str(ROOT / "license.dat"),
                  str(Path.cwd() / "license.dat"), hx.HxVersion.license_path]
    for c in candidates:
        if c and Path(c).is_file():
            hx.HxVersion.set_license_path(c)
            return c
    raise SystemExit(
        "Hexaly licence not found. Download license.dat from your hexaly.com account and put it in "
        f"{ROOT} (or pass --license PATH, or set HX_LICENSE_PATH).")


def solve(
    inst: Instance,
    prm: Params,
    n_stops: int,
    stops_mode: str = "eq",
    time_limit_s: int | None = 60,
    random_seed: int = 0,
    threads: int | None = None,
    use_greedy: bool = False,
    verbose: bool = True,
    iteration_limit: int | None = None,
    stall_s: float | None = None,
    target: float | None = None,
) -> dict:
    t_start = time.perf_counter()
    res = {"solver": "Hexaly", "n_stops": n_stops, "stops_mode": stops_mode, "random_seed": random_seed,
           "greedy_start": use_greedy}
    S = position_sets(inst, n_stops, stops_mode)
    keys = sorted(set().union(*S[1:]))
    if not S[1] or (stops_mode == "eq" and not S[n_stops]):
        return {**res, "status": "INFEASIBLE", "build_time_s": 0.0, "solve_time_s": 0.0}
    idx = {k: i for i, k in enumerate(keys)}
    n = len(keys)

    # --- data as dense arrays (Hexaly "at" operator needs arrays)
    gate_t = [MISSING_ARC_T] * n
    gate_d = [MISSING_ARC_D] * n
    gate_miss = [1] * n
    T = [[MISSING_ARC_T] * n for _ in range(n)]
    D = [[MISSING_ARC_D] * n for _ in range(n)]
    M = [[1] * n for _ in range(n)]
    for (k, l), (t, d) in inst.arcs.items():
        if l not in idx:
            continue
        if k == GATE_ID:
            gate_t[idx[l]], gate_d[idx[l]], gate_miss[idx[l]] = t, d, 0
        elif k in idx:
            a, b = idx[k], idx[l]
            T[a][b], D[a][b], M[a][b] = t, d, 0
    C = [[inst.poi_cost[k].get(j, MISSING_POI) for j in inst.pois] for k in keys]
    demand = [inst.demand[j] for j in inst.pois]

    with hx.HexalyOptimizer() as opt:
        m = opt.model
        gT, gD, gM = m.array(gate_t), m.array(gate_d), m.array(gate_miss)
        aT, aD, aM, aC = m.array(T), m.array(D), m.array(M), m.array(C)

        route = m.list(n)
        c = m.count(route)
        if stops_mode == "eq":
            m.constraint(c == n_stops)
        else:
            m.constraint(c >= 1)
            m.constraint(c <= n_stops)

        missing = gM[route[0]] + m.sum(m.range(1, c), m.lambda_function(
            lambda i: m.at(aM, route[i - 1], route[i])))
        arrival = m.array(m.range(0, c), m.lambda_function(
            lambda i, prev: m.iif(i == 0, gT[route[0]], prev + m.at(aT, route[i - 1], route[i]))), 0)
        t_city = arrival[c - 1]
        length = gD[route[0]] + m.sum(m.range(1, c), m.lambda_function(
            lambda i: m.at(aD, route[i - 1], route[i])))
        buses = m.ceil((prm.t_outside_h + t_city / 60 + prm.t_rest_h) * prm.frequency - 1e-9)

        def best_time(j: int):  # factory: Hexaly counts lambda arguments, so no default args
            return m.min(m.range(0, c), m.lambda_function(
                lambda p: prm.w_iu * arrival[p] + m.at(aC, route[p], j)))

        user_min = m.sum([demand[j] * best_time(j) for j in range(len(demand))])
        c_op = prm.ct * prm.t0 * buses + prm.cd * length
        c_user = prm.vot / 60 * user_min
        cost = prm.rho_operator * c_op + prm.rho_user * c_user

        m.minimize(missing)
        m.minimize(cost)
        m.close()

        if use_greedy:
            g = greedy_route(inst, prm, S, n_stops, stops_mode)
            if g:
                route.value.clear()
                for k in g:
                    route.value.add(idx[k])
                res["hint_objective"] = evaluate_route(inst, prm, g)["objective"]

        if not time_limit_s and not iteration_limit:
            raise ValueError("set a time limit or an iteration limit (stall rules alone may never stop)")
        if time_limit_s:
            opt.param.time_limit = int(math.ceil(time_limit_s))
        if iteration_limit:
            opt.param.iteration_limit = iteration_limit
        if target is not None:
            opt.param.set_objective_threshold(0, 0)       # no missing arcs ...
            opt.param.set_objective_threshold(1, target)  # ... and cost <= target
        opt.param.seed = random_seed
        if threads:
            opt.param.nb_threads = threads
        opt.param.verbosity = 1 if verbose else 0
        opt.param.time_between_ticks = 1

        # best feasible cost over time (1 s resolution)
        trace = {"best": math.inf, "first": None, "best_t": None, "stopped": None}
        t_solve = [0.0]

        def on_tick(o, cb_type):
            # before the first feasible solution some expressions (e.g. route[0]) are undefined
            if o.solution.status not in (hx.HxSolutionStatus.FEASIBLE, hx.HxSolutionStatus.OPTIMAL):
                return
            try:
                if missing.value != 0:
                    return
                v = cost.value
            except hx.HxError:
                return
            t = time.perf_counter() - t_solve[0]
            if trace["first"] is None:
                trace["first"] = t
            if v < trace["best"] - 1e-6:
                trace["best"], trace["best_t"] = v, t
            if stall_s and t - trace["best_t"] >= stall_s:
                trace["stopped"] = "stall_time"
                o.stop()

        opt.add_callback(hx.HxCallbackType.TIME_TICKED, on_tick)
        build_s = time.perf_counter() - t_start
        if verbose:
            print(f"[Hexaly N={n_stops}] {n} stops (built in {build_s:.1f}s)", flush=True)
        t_solve[0] = time.perf_counter()
        opt.solve()
        solve_s = time.perf_counter() - t_solve[0]

        sol = opt.solution
        status = str(sol.status).split(".")[-1]
        n_iter = opt.statistics.nb_iterations
        if trace["stopped"]:
            stop_reason = trace["stopped"]
        elif target is not None and trace["best"] <= target + 1e-6:
            stop_reason = "target"
        elif iteration_limit and n_iter >= iteration_limit:
            stop_reason = "iterations"
        elif str(sol.status).endswith("OPTIMAL"):
            stop_reason = "optimal"
        else:
            stop_reason = "time_limit"
        res["stop_reason"] = stop_reason
        if sol.status not in (hx.HxSolutionStatus.FEASIBLE, hx.HxSolutionStatus.OPTIMAL):
            res.update({"build_time_s": build_s, "solve_time_s": solve_s, "hexaly_status": status,
                        "n_iterations": n_iter, "status": "INFEASIBLE"})
            return res
        route_idx = [int(v) for v in route.value]
        n_missing = int(missing.value)
        bound = sol.get_objective_bound(1)
        hx_cost = cost.value
        if target is not None and n_missing == 0 and hx_cost <= target + 1e-6 and res["stop_reason"] == "time_limit":
            res["stop_reason"] = "target"

    res.update({"build_time_s": build_s, "solve_time_s": solve_s, "hexaly_status": status, "n_iterations": n_iter})
    route_keys = [keys[i] for i in route_idx]
    if n_missing > 0 or not route_keys:
        res.update({"status": "INFEASIBLE", "route": route_keys})
        return res
    ev = evaluate_route(inst, prm, route_keys)
    res.update({
        "status": "OPTIMAL" if status == "OPTIMAL" else "FEASIBLE",
        "hexaly_objective": hx_cost,
        "best_bound": bound,
        "gap": abs(hx_cost - bound) / max(abs(hx_cost), 1e-9),
        "time_to_first_s": trace["first"],
        "time_to_best_s": trace["best_t"],
        "route_names": [inst.names.get(k, k) for k in route_keys],
        **ev,
    })
    if abs(ev["objective"] - hx_cost) > 1e-4 * abs(ev["objective"]):
        print(f"  WARNING: Hexaly objective {hx_cost:.4f} != evaluated objective {ev['objective']:.4f}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stops", type=int, nargs="+", default=[2])
    ap.add_argument("--stops-mode", choices=["eq", "max"], default="eq")
    ap.add_argument("--time-limit", type=int, default=60,
                    help="max seconds per run, integer (0 = no time cap; then set --iterations)")
    ap.add_argument("--iterations", type=int, default=None, help="max Hexaly iterations (moves tried)")
    ap.add_argument("--stall-time", type=float, default=None,
                    help="stop after this many seconds without improvement (1 s resolution)")
    ap.add_argument("--target", type=float, default=None, help="stop once the objective is <= this value")
    ap.add_argument("--rs", type=int, nargs="+", default=[0], help="Hexaly random seed(s)")
    ap.add_argument("--threads", type=int, default=None, help="Hexaly threads (default: Hexaly's choice)")
    ap.add_argument("--greedy", action="store_true", help="start from the greedy route")
    ap.add_argument("--license", default=None, help="path to license.dat")
    ap.add_argument("--check", action="store_true", help="compare with brute-force enumeration (small instances)")
    ap.add_argument("--quiet", action="store_true", help="hide the Hexaly log")
    ap.add_argument("--out-dir", default=str(ROOT / "results/hexaly"))
    ap.add_argument("--sample-stops", type=int, default=None, help="random subset of candidate stops (testing)")
    ap.add_argument("--seed", type=int, default=0, help="seed of the random stop subset")
    prm_default = Params()
    for name, val in vars(prm_default).items():
        ap.add_argument(f"--{name.replace('_', '-')}", type=type(val), default=val)
    args = ap.parse_args()

    lic = set_license(args.license)
    print(f"Hexaly {hx.HxVersion.version}, licence: {lic}", flush=True)
    prm = Params(**{k: getattr(args, k) for k in vars(prm_default)})
    print("Loading data...", flush=True)
    inst = load_instance(sample_stops=args.sample_stops, seed=args.seed)
    print(f"{len(inst.stops)} candidate stops, {len(inst.arcs)} arcs, {len(inst.pois)} POIs", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = f"_S{args.sample_stops}_seed{args.seed}" if args.sample_stops else ""
    summary = []
    for n in args.n_stops:
        bf = None
        if args.check:
            from b0_check_exact import brute_force
            bf = brute_force(inst, prm, n, args.stops_mode)[0]
        for rs in args.rs:
            res = solve(inst, prm, n, args.stops_mode, args.time_limit or None, rs, args.threads, args.greedy,
                        not args.quiet, args.iterations, args.stall_time, args.target)
            res.update({"sample_stops": args.sample_stops, "seed": args.seed, "time_limit_s": args.time_limit or None,
                        "max_iterations": args.iterations, "stall_time_s": args.stall_time, "target": args.target,
                        "threads": args.threads})
            (out_dir / f"Hexaly_N{n}_{args.stops_mode}{tag}_rs{rs}.json").write_text(
                json.dumps({"params": vars(prm), **res}, indent=2, ensure_ascii=False))
            summary.append(res)
            if "objective" in res:
                print(f"  {res['status']} obj={res['objective']:.2f} (bound {res['best_bound']:.2f}) "
                      f"C_O={res['operator_cost']:.2f} C_U={res['user_cost']:.2f} buses={res['buses']} "
                      f"T_C={res['city_time_min']:.1f}min best found at ~{res['time_to_best_s'] or 0:.0f}s "
                      f"/ {res['solve_time_s']:.1f}s, {res['n_iterations']} iterations, "
                      f"stopped by {res['stop_reason']}")
                print("  route: Gate -> " + " -> ".join(res["route_names"]))
                if bf is not None:
                    gap = 100 * (res["objective"] / bf - 1)
                    print(f"  brute force optimum {bf:.2f}: " + ("OK" if abs(gap) <= 0.01 else f"gap {gap:.3f}%"))
            else:
                print(f"  {res['status']} after {res['solve_time_s']:.1f}s")

    cols = ["solver", "sample_stops", "seed", "random_seed", "n_stops", "stops_mode", "status", "objective",
            "best_bound", "gap", "operator_cost", "user_cost", "buses", "city_time_min", "length_km",
            "avg_user_time_min", "time_to_first_s", "time_to_best_s", "solve_time_s", "stop_reason", "n_iterations",
            "time_limit_s", "max_iterations", "stall_time_s", "target", "threads",
            "greedy_start", "route"]
    df = pd.DataFrame([{c: r.get(c) for c in cols} for r in summary])
    df["route"] = df["route"].apply(lambda r: " -> ".join(r) if isinstance(r, list) else r)
    df.to_csv(out_dir / "summary.csv", index=False, mode="a", header=not (out_dir / "summary.csv").exists())
    print(df.drop(columns="route").to_string(index=False))


if __name__ == "__main__":
    main()
