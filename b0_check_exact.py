"""Validate the exact MILP (b1_exact_ortools.py) on random subsets of candidate stops.

For each (subset size, seed, N), every route Gate -> k1 -> ... -> kN is enumerated by
brute force and evaluated with `evaluate_route` (independent of the MILP). The optimum
must match both the MILP objective and the re-evaluated cost of the MILP route.

python b0_check_exact.py                                   # default grid of checks
python b0_check_exact.py --sizes 40 60 --seeds 0 1 2 --n-stops 1 2 3 --solvers SCIP HIGHS CBC
python b0_check_exact.py --solvers CPSAT HIGHS                # CP-SAT model of b2_cpsat.py
"""

from __future__ import annotations

import argparse
import time

from b1_exact_ortools import build_and_solve
from b2_cpsat import build_and_solve as cpsat_solve
from bndp_common import GATE_ID, Instance, Params, evaluate_route, load_instance, subset_instance


def brute_force(inst: Instance, prm: Params, n_stops: int, stops_mode: str) -> tuple[float, list[str], int]:
    out = inst.out_arcs()
    best = (float("inf"), [], 0)
    count = 0

    def dfs(route: list[str], visited: set[str]) -> None:
        nonlocal best, count
        if route and (len(route) == n_stops or stops_mode == "max"):
            count += 1
            obj = evaluate_route(inst, prm, route)["objective"]
            if obj < best[0]:
                best = (obj, list(route), 0)
        if len(route) == n_stops:
            return
        for l, _, _ in out.get(route[-1] if route else GATE_ID, []):
            if l not in visited:
                route.append(l)
                visited.add(l)
                dfs(route, visited)
                visited.discard(l)
                route.pop()

    dfs([], set())
    return best[0], best[1], count


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=[30, 60])
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2])
    ap.add_argument("--n-stops", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--solvers", nargs="+", default=["SCIP", "HIGHS"])
    ap.add_argument("--stops-mode", choices=["eq", "max"], default="eq")
    args = ap.parse_args()

    prm = Params()
    full = load_instance()
    failures, rows = 0, []
    for size in args.sizes:
        for seed in args.seeds:
            inst = subset_instance(full, size, seed)
            for n in args.n_stops:
                t0 = time.perf_counter()
                bf_obj, bf_route, n_routes = brute_force(inst, prm, n, args.stops_mode)
                bf_t = time.perf_counter() - t0
                line = f"S={size:<4} seed={seed} N={n}  brute force: {n_routes:>7} routes, obj={bf_obj:10.3f} ({bf_t:5.1f}s)"
                for s in args.solvers:
                    if s.upper().startswith("CPSAT"):
                        res = cpsat_solve(inst, prm, n, args.stops_mode, time_limit_s=120,
                                          lns_only=s.upper() == "CPSAT-LNS", verbose=False)
                    else:
                        res = build_and_solve(inst, prm, n, s, args.stops_mode, rel_gap=1e-6)
                    obj = res.get("objective", float("inf"))
                    mip = res.get("mip_objective", obj)
                    tol = 1e-4 * max(1, abs(bf_obj))
                    ok = (res["status"] == "INFEASIBLE" and n_routes == 0) or (
                        abs(obj - bf_obj) <= tol and abs(mip - bf_obj) <= tol
                    )
                    failures += not ok
                    line += f" | {s}: {res['status'][:4]} {obj:10.3f} {res['solve_time_s']:5.1f}s {'OK' if ok else 'MISMATCH'}"
                    if not ok:
                        line += f"\n    brute={bf_route} milp={res.get('route')}"
                rows.append(line)
                print(line, flush=True)

    print("\n".join(["", "=" * 40, *rows, "=" * 40]))
    print("ALL CHECKS PASSED" if failures == 0 else f"{failures} MISMATCHES")
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
