"""Exact MILP for the interurban Bus Network Design Problem (Roca-Riu et al., 2012) with OR-Tools.

Single interurban route that enters Barcelona at the city gate and visits exactly
(or at most) N candidate stops. A single path has exactly one terminal, so Eq. (11)
|P| <= p holds for any p >= 1.

Bilevel -> single level
-----------------------
The lower level is an All-or-Nothing assignment: every POI demand alights at the
route stop minimising its own perceived time. Since the upper level minimises
rho_O*C_O + rho_U*C_U with rho_U > 0, C_U is separable by POI and C_O does not
depend on the assignment, minimising over the assignment variables jointly gives
exactly the users' best responses. Hence the single-level MILP below is exact.

Position-based (layered) formulation
------------------------------------
Positions p = 1..N along the route; S_p = stops that can be the p-th stop
(reachable from the gate in exactly p arcs and, for exactly-N routes, able to
continue for N-p more arcs). A_p = arcs from position p-1 to p (p=1 starts at the gate).

x_klp in {0,1}  arc (k,l) used between positions p-1 and p    -- the only binaries
u_kp  in [0,1]  stop k at position p  (= inflow of x at (k,p))
T_p   >= 0      in-city travel time of the arc entering position p (min)
z_jkp in [0,1]  POI j alights at stop k at position p (lower level)
h_jp  in {0,1}  = sum_k z_jkp: POI j alights at position p
G_jp  = sum_{p'>=p} h_jp' (POI j still on board on arc p)
w_jp  >= 0      = T_p * G_jp (McCormick, exact because G is binary); G_j1 = 1
b     in Z+     number of buses

min rho_O * (ct*t0*b + cd*sum d_kl x_klp)                                 (Eqs. 1-5)
  + rho_U * VoT/60 * sum_j d_j (w_iu*(T_1 + sum_{p>=2} w_jp) + sum c_kj z_jkp)   (Eq. 6)
s.t. sum x_{gate,l,1} = 1                                                 (Eqs. 8-10)
     sum_l x_{k,l,p+1} = u_kp  (exactly N)  or  <= u_kp  (at most N)
     sum_p u_kp <= 1                                     (each stop visited once)
     T_p = sum t_kl x_klp;   b >= f (T_A + T_rest + sum_p T_p / 60)        (Eq. 3)
     sum_kp z_jkp = 1;  z_jkp <= u_kp
     w_jp >= T_p - Tmax_p (1 - G_jp);  w_jp >= Tmin_p G_jp

Positions make the route a path in a layered DAG, so no subtour / arrival-time
big-M constraints are needed, and arrival times are linear in x. h must be binary:
with a fractional split between positions the McCormick term would under-count the
in-vehicle time. Given the route and h, each position holds a single stop, so z can
stay continuous and the model is exact.

Usage
-----
python b1_exact_ortools.py --n-stops 3 5 8 --time-limit 600          # HiGHS by default
python b1_exact_ortools.py --n-stops 6 --solver HIGHS SCIP --stops-mode max
python b1_exact_ortools.py --n-stops 2 3 4 --solver SCIP HIGHS --sample-stops 80 --seed 1   # quick test
Validation against brute-force enumeration on random subsets: python b0_check_exact.py
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from ortools.linear_solver import pywraplp

from bndp_common import GATE_ID, ROOT, Instance, Params, evaluate_route, load_instance

STATUS = {
    pywraplp.Solver.OPTIMAL: "OPTIMAL",
    pywraplp.Solver.FEASIBLE: "FEASIBLE",
    pywraplp.Solver.INFEASIBLE: "INFEASIBLE",
    pywraplp.Solver.UNBOUNDED: "UNBOUNDED",
    pywraplp.Solver.ABNORMAL: "ABNORMAL",
    pywraplp.Solver.NOT_SOLVED: "NOT_SOLVED",
}


def position_sets(inst: Instance, n_stops: int, stops_mode: str) -> list[set[str]]:
    """S_p for p = 1..N (index 0 unused): stops that can occupy position p."""
    out = inst.out_arcs()
    fwd: list[set[str]] = [set(), {l for l, _, _ in out.get(GATE_ID, [])}]
    for _ in range(2, n_stops + 1):
        fwd.append({l for k in fwd[-1] for l, _, _ in out.get(k, [])})
    if stops_mode == "max":
        return fwd
    pos = [set() for _ in range(n_stops + 1)]
    pos[n_stops] = fwd[n_stops]
    for p in range(n_stops - 1, 0, -1):
        pos[p] = {k for k in fwd[p] if any(l in pos[p + 1] for l, _, _ in out.get(k, []))}
    return pos


def build_and_solve(
    inst: Instance,
    prm: Params,
    n_stops: int,
    solver_name: str = "SCIP",
    stops_mode: str = "eq",
    time_limit_s: float | None = None,
    rel_gap: float | None = None,
    threads: int | None = None,
    verbose: bool = False,
) -> dict:
    t_start = time.perf_counter()
    solver = pywraplp.Solver.CreateSolver(solver_name)
    if solver is None:
        raise RuntimeError(f"Solver {solver_name} not available in this OR-Tools build")
    res = {"solver": solver_name, "n_stops": n_stops, "stops_mode": stops_mode}

    # --- Layered network: S_p and arcs between consecutive positions
    S = position_sets(inst, n_stops, stops_mode)
    P = [p for p in range(1, n_stops + 1) if S[p]]
    if not S[1] or (stops_mode == "eq" and len(P) < n_stops):
        return {**res, "status": "INFEASIBLE", "build_time_s": 0.0, "solve_time_s": 0.0,
                "n_vars": 0, "n_constraints": 0}
    out = inst.out_arcs()
    layer_arcs: dict[int, dict[tuple[str, str], tuple[float, float]]] = {}
    for p in P:
        tails = [GATE_ID] if p == 1 else sorted(S[p - 1])
        layer_arcs[p] = {(k, l): (t, d) for k in tails for l, t, d in out.get(k, []) if l in S[p]}

    # --- Route variables
    x = {p: {a: solver.BoolVar(f"x[{a[0]}->{a[1]},{p}]") for a in layer_arcs[p]} for p in P}
    u = {(k, p): solver.NumVar(0, 1, f"u[{k},{p}]") for p in P for k in S[p]}
    in_x = {kp: [] for kp in u}
    out_x = {kp: [] for kp in u}
    for p in P:
        for (k, l), var in x[p].items():
            in_x[(l, p)].append(var)
            if p > 1:
                out_x[(k, p - 1)].append(var)

    solver.Add(sum(x[1].values()) == 1)
    for (k, p), var in u.items():
        solver.Add(var == sum(in_x[(k, p)]))
        if p + 1 in x:
            if stops_mode == "eq":
                solver.Add(sum(out_x[(k, p)]) == var)
            else:
                solver.Add(sum(out_x[(k, p)]) <= var)
    positions_of: dict[str, list] = {}
    for (k, p), var in u.items():
        positions_of.setdefault(k, []).append(var)
    for k, vs in positions_of.items():
        if len(vs) > 1:
            solver.Add(sum(vs) <= 1)

    T = {p: solver.NumVar(0, max(t for t, _ in layer_arcs[p].values()), f"T[{p}]") for p in P}
    t_min = {p: min(t for t, _ in layer_arcs[p].values()) for p in P}
    t_max = {p: max(t for t, _ in layer_arcs[p].values()) for p in P}
    for p in P:
        solver.Add(T[p] == sum(t * x[p][a] for a, (t, _) in layer_arcs[p].items()))
    length = sum(d * x[p][a] for p in P for a, (_, d) in layer_arcs[p].items())

    b = solver.IntVar(0, solver.infinity(), "buses")
    solver.Add(b >= prm.frequency * (prm.t_outside_h + prm.t_rest_h + sum(T.values()) / 60))  # Eq. 3

    # --- Lower level: POI j alights at stop k in position p
    user_terms = []
    for j in inst.pois:
        h = {}
        egress = []
        for p in P:
            zs = []
            for k in S[p]:
                c = inst.poi_cost[k].get(j)
                if c is None:
                    continue
                z = solver.NumVar(0, 1, f"z[{j},{k},{p}]")
                solver.Add(z <= u[(k, p)])
                zs.append(z)
                egress.append(c * z)
            h[p] = solver.BoolVar(f"h[{j},{p}]")
            solver.Add(h[p] == sum(zs))
        solver.Add(sum(h.values()) == 1)

        ride = T[P[0]]
        for idx, p in enumerate(P[1:], start=1):
            G = sum(h[q] for q in P[idx:])
            w = solver.NumVar(0, t_max[p], f"w[{j},{p}]")
            solver.Add(w >= T[p] - t_max[p] * (1 - G))
            solver.Add(w >= t_min[p] * G)
            ride = ride + w
        user_terms.append(inst.demand[j] * (prm.w_iu * ride + sum(egress)))

    c_op = prm.ct * prm.t0 * b + prm.cd * length
    c_user = prm.vot / 60 * sum(user_terms)
    solver.Minimize(prm.rho_operator * c_op + prm.rho_user * c_user)

    # --- Solve
    if verbose:
        solver.EnableOutput()
    if time_limit_s:
        solver.SetTimeLimit(int(time_limit_s * 1000))
    if threads:
        solver.SetNumThreads(threads)
    sp = pywraplp.MPSolverParameters()
    if rel_gap is not None:
        sp.SetDoubleParam(pywraplp.MPSolverParameters.RELATIVE_MIP_GAP, rel_gap)

    build_s = time.perf_counter() - t_start
    n_arcs = sum(len(a) for a in layer_arcs.values())
    print(
        f"[{solver_name} N={n_stops}] {solver.NumVariables()} vars, {solver.NumConstraints()} constraints, "
        f"{len(positions_of)} stops, {n_arcs} layered arcs (built in {build_s:.1f}s)",
        flush=True,
    )
    status = solver.Solve(sp)
    solve_s = time.perf_counter() - t_start - build_s

    res.update({
        "status": STATUS.get(status, str(status)),
        "build_time_s": build_s,
        "solve_time_s": solve_s,
        "n_vars": solver.NumVariables(),
        "n_constraints": solver.NumConstraints(),
    })
    if status not in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE):
        return res

    # --- Extract route: the arc used at each position
    route = []
    for p in P:
        used = [l for (_, l), var in x[p].items() if var.solution_value() > 0.5]
        if not used:
            break
        route.append(used[0])

    best_bound = solver.Objective().BestBound()
    obj = solver.Objective().Value()
    ev = evaluate_route(inst, prm, route)
    res.update(
        {
            "mip_objective": obj,
            "best_bound": best_bound,
            "gap": abs(obj - best_bound) / max(abs(obj), 1e-9),
            "route_names": [inst.names.get(k, k) for k in route],
            **ev,
        }
    )
    if abs(ev["objective"] - obj) > 1e-4 * max(1.0, abs(obj)):
        print(f"  WARNING: MIP objective {obj:.4f} != evaluated objective {ev['objective']:.4f}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-stops", type=int, nargs="+", default=[2], help="number of stops (one solve per value)")
    ap.add_argument("--stops-mode", choices=["eq", "max"], default="eq", help="exactly N stops or at most N")
    ap.add_argument("--solver", nargs="+", default=["HIGHS"], help="HIGHS (fastest here), SCIP, CBC, ...")
    ap.add_argument("--time-limit", type=float, default=3600*3, help="seconds per solve")
    ap.add_argument("--gap", type=float, default=0.01, help="relative MIP gap")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--verbose", action="store_true", help="show solver log")
    ap.add_argument("--out-dir", default=str(ROOT / "results/exact"))
    ap.add_argument("--sample-stops", type=int, default=None,
                    help="use a random subset of this many candidate stops (for testing)")
    ap.add_argument("--seed", type=int, default=0, help="seed of the random stop subset")
    # model parameters
    prm_default = Params()
    for name, val in vars(prm_default).items():
        ap.add_argument(f"--{name.replace('_', '-')}", type=type(val), default=val)
    args = ap.parse_args()

    prm = Params(**{k: getattr(args, k) for k in vars(prm_default)})
    print("Loading data...", flush=True)
    inst = load_instance(sample_stops=args.sample_stops, seed=args.seed)
    print(f"{len(inst.stops)} candidate stops, {len(inst.arcs)} arcs, {len(inst.pois)} POIs, "
          f"{sum(inst.demand.values()):.0f} trips", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    tag = f"_S{args.sample_stops}_seed{args.seed}" if args.sample_stops else ""
    for solver_name in args.solver:
        for n in args.n_stops:
            res = build_and_solve(
                inst, prm, n, solver_name, args.stops_mode, args.time_limit, args.gap, args.threads, args.verbose
            )
            res["sample_stops"], res["seed"] = args.sample_stops, args.seed
            (out_dir / f"{solver_name}_N{n}_{args.stops_mode}{tag}.json").write_text(
                json.dumps({"params": vars(prm), **res}, indent=2, ensure_ascii=False)
            )
            summary.append(res)
            if "objective" in res:
                print(
                    f"  {res['status']} obj={res['objective']:.2f} (bound {res['best_bound']:.2f}, "
                    f"gap {100 * res['gap']:.2f}%) C_O={res['operator_cost']:.2f} C_U={res['user_cost']:.2f} "
                    f"buses={res['buses']} T_C={res['city_time_min']:.1f}min L={res['length_km']:.2f}km "
                    f"time={res['solve_time_s']:.1f}s"
                )
                print("  route: Gate -> " + " -> ".join(res["route_names"]))
            else:
                print(f"  {res['status']} after {res['solve_time_s']:.1f}s")

    cols = ["solver", "sample_stops", "seed", "n_stops", "status", "objective", "best_bound", "gap", "operator_cost", "user_cost",
            "buses", "city_time_min", "length_km", "avg_user_time_min", "solve_time_s", "route"]
    import pandas as pd

    df = pd.DataFrame([{c: r.get(c) for c in cols} for r in summary])
    df["route"] = df["route"].apply(lambda r: " -> ".join(r) if isinstance(r, list) else r)
    df.to_csv(out_dir / "summary.csv", index=False, mode="a", header=not (out_dir / "summary.csv").exists())
    print(df.drop(columns="route").to_string(index=False))


if __name__ == "__main__":
    main()
