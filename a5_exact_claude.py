"""
Interurban Bus Network Design Problem (BNDP) -- Roca-Riu, Estrada & Trapote (2012)
Single route, OR-Tools MathOpt.  Equation numbers refer to the paper.

Written from the PuLP model (interurban_bnd_pulp.py), with ONE change in the user cost:

  Eq. (6) charges every user t_ik,IU, the time with interurban service between the city gate
  and the stop k where the user gets off. A bus only reaches k by driving the route, so that
  time is the sum of the stop-to-stop times t_kl,IU along the route from the gate to k (the
  same arc times that Eq. (4) adds up for the whole route). The cost of a user who gets off at
  stop k and goes to destination j is therefore

        VoT * d_j * ( T_k + t_kj )          T_k = arrival time at k along the route

  where t_kj is the local chain (transfer + waiting + urban ride + egress, already weighted).
  The PuLP model used the direct gate->k edge instead, which ignores the route.

How T_k is written in the MILP (T_k is only needed for the stops that are on the route):
  * T_gate = 0 and, for every used arc (k, l):  T_l >= T_k + t_kl          (arrival times)
  * w_kj >= T_k - T_cap*(1 - z_kj)       so w_kj = T_k when z_kj = 1 and 0 otherwise
    The cost uses w_kj, so the solver keeps T_k at the smallest value the arcs allow, which is
    the real cumulative time.  z must stay binary: this link is only exact for z in {0,1}.
  * Every arc takes > 0 min, so a cycle would need T to grow around it: loops are impossible.
    The T variables therefore replace the MTZ constraints of the PuLP model.

NOT in the paper (marked in the code): T_cap, a cap on the minutes spent inside the city, and the
optional number-of-stops constraint.
"""

from __future__ import annotations

import argparse
from datetime import timedelta
import json
import math
from pathlib import Path

import pandas as pd
from ortools.math_opt.python import mathopt


# ---------------------------------------------------------------------------------------------
# Evaluation of one given route (plain Python). It is used to report the solution and to check
# the MILP: at the optimum, the MILP objective must equal this value.
# ---------------------------------------------------------------------------------------------
def evaluate_route(route, gate, J, demand, t_min, dist, *, freq, t_outside_h, t_rest_h,
                   c_d, c_t, t0, vot, rho_o=0.5, rho_u=0.5):
    """route = [gate, k1, k2, ...]. Users choose their best drop-off stop on the route (AoN)."""
    cum, elapsed, length_km = {gate: 0.0}, 0.0, 0.0
    for a, b in zip(route, route[1:]):
        elapsed += t_min[(a, b)]          # Eq. (4)
        length_km += dist[(a, b)]         # Eq. (5)
        cum[b] = elapsed                  # T_b: arrival time at b
    time_in_city_h = elapsed / 60.0
    buses = math.ceil((t_outside_h + time_in_city_h + t_rest_h) * freq - 1e-9)       # Eq. (3)
    operator_cost = buses * c_t * t0 + c_d * length_km                                # Eq. (1)

    transfer_stop, user_hours = {}, 0.0
    for j in J:
        candidates = [k for k in route if (k, j) in t_min]
        k_best = min(candidates, key=lambda k: cum[k] + t_min[(k, j)])
        transfer_stop[j] = k_best
        user_hours += demand[j] * (cum[k_best] + t_min[(k_best, j)]) / 60.0
    user_cost = vot * user_hours                                                      # Eq. (6)
    return {
        "cum": cum,
        "transfer_stop": transfer_stop,
        "buses": buses,
        "time_in_city_h": time_in_city_h,
        "length_km": length_km,
        "operator_cost": operator_cost,
        "user_cost": user_cost,
        "objective": rho_o * operator_cost + rho_u * user_cost,                       # Eq. (7)
    }


def build_model(K, J, gate, demand, t_min, dist, *, freq, t_outside_h, t_rest_h, c_d, c_t, t0,
                vot, rho_o=0.5, rho_u=0.5, p_max=1, n_stops=None, stops_mode="max",
                t_max_min=85.0):
    """
    K        : candidate stops, including the gate
    J        : destinations
    demand   : dict destination -> d_j
    t_min    : dict (a, b) -> minutes. (stop, stop) and (gate, stop): in-bus time t_kl,IU of the arc.
                              (stop, destination): local chain time, already weighted.
    dist     : dict (a, b) -> km, for the (stop, stop) and (gate, stop) arcs
    t_max_min: cap on the minutes spent inside the city, one way        [NOT in the paper]
               When n_stops is given, a tighter valid cap is derived from the arc times and the
               smaller of the two is used. None = rely only on n_stops.
    n_stops  : number of stops after the gate (None = free length, as in the paper)  [NOT in the paper]
    stops_mode: "exact" or "max"
    """
    kid = {k: i for i, k in enumerate(K)}
    jid = {j: i for i, j in enumerate(J)}
    print(f'len candidate stops K = {len(K)}, destinations J = {len(J)}')

    # arcs into the gate are never used by a route, so they are not created
    arcs = [(k, l) for k in K for l in K
            if k != l and l != gate and (k, l) in t_min and (k, l) in dist]

    print(f'len arcs = {len(arcs)}')
    if not any(a[0] == gate for a in arcs):
        raise ValueError("No arcs leave the city gate")
    if any(t_min[a] <= 0 for a in arcs):
        raise ValueError("Every stop->stop arc must take > 0 min (needed to exclude loops)")
    pairs = [(k, j) for k in K for j in J if (k, j) in t_min]
    print(f"Feasible stop-destination pairs: {len(pairs)}")
    uncovered = set(J) - {j for _, j in pairs}
    if uncovered:
        raise ValueError(f"Destinations without any reachable transfer stop: {sorted(uncovered)}")

    # valid upper bound of the time in the city; it is also the big-M of the T constraints
    gate_max = max(t_min[a] for a in arcs if a[0] == gate)
    inner = [t_min[a] for a in arcs if a[0] != gate]
    t_cap = t_max_min
    if n_stops is not None:
        bound = gate_max + (n_stops - 1) * (max(inner) if inner else 0.0)
        t_cap = bound if t_cap is None else min(t_cap, bound)
    if t_cap is None:
        raise ValueError("Pass t_max_min or n_stops: the arrival times need an upper bound")

    out_arcs = {k: [] for k in K}
    in_arcs = {k: [] for k in K}
    for a in arcs:
        out_arcs[a[0]].append(a)
        in_arcs[a[1]].append(a)

    model = mathopt.Model(name="Interurban_BNDP_cumulative_time")

    # ---- decision variables
    x = {a: model.add_binary_variable(name=f"X_{kid[a[0]]}_{kid[a[1]]}") for a in arcs}   # (15)
    term = {k: model.add_binary_variable(name=f"Term_{kid[k]}") for k in K if k != gate}   # k in P
    z = {p: model.add_binary_variable(name=f"Z_{kid[p[0]]}_{jid[p[1]]}") for p in pairs}   # transfer stop
    buses = model.add_integer_variable(lb=0, name="Buses")
    # arrival time (min) at each stop, measured from the gate
    T = {k: model.add_variable(lb=0.0, name=f"T_{kid[k]}") for k in K if k != gate}
    #T = {k: model.add_variable(lb=0.0, ub=t_cap, name=f"T_{kid[k]}") for k in K if k != gate}

    # w_kj = in-bus time of the users that get off at k and go to j
    w = {(k, j): model.add_variable(lb=0.0, name=f"W_{kid[k]}_{jid[j]}")
         for k, j in pairs if k != gate}

    # w = {(k, j): model.add_variable(lb=0.0, ub=t_cap, name=f"W_{kid[k]}_{jid[j]}")
    #          for k, j in pairs if k != gate}
    # on_route_k = 1 if the route stops at k. It equals the sum of the (at most one) arcs entering k,
    # so it is integral without being declared binary.
    on_route = {k: model.add_variable(lb=0.0, ub=1.0, name=f"On_{kid[k]}") for k in K if k != gate}

    # ---- operator cost: Eqs. (1)-(5)
    time_in_city = mathopt.fast_sum((t_min[a] / 60.0) * x[a] for a in arcs)     # (4), hours
    length_km = mathopt.fast_sum(dist[a] * x[a] for a in arcs)                  # (5)
    operator_cost = buses * c_t * t0 + c_d * length_km                          # (1)
    # (3): buses >= cycle time * frequency. `buses` has a positive cost, so at the optimum it is
    # the smallest integer that satisfies this, i.e. the ceiling.
    model.add_linear_constraint(buses >= (t_outside_h + time_in_city + t_rest_h) * freq,
                                name="Eq3_buses")

    # ---- user cost: Eq. (6) with the in-bus time accumulated along the route
    user_cost = vot * (
        mathopt.fast_sum(demand[j] * (t_min[(k, j)] / 60.0) * z[(k, j)] for k, j in pairs)
        + mathopt.fast_sum(demand[j] / 60.0 * w[(k, j)] for k, j in w)
    )

    # ---- objective: Eq. (7)
    model.minimize(rho_o * operator_cost + rho_u * user_cost)

    # ---- route structure: Eqs. (8)-(11)
    for k in K:
        out_k = mathopt.fast_sum(x[a] for a in out_arcs[k])
        in_k = mathopt.fast_sum(x[a] for a in in_arcs[k])
        if k == gate:
            model.add_linear_constraint(out_k - in_k == 1, name="Eq9_route_starts_at_gate")
        else:
            # term_k = 0: out - in = 0 (8);   term_k = 1: the route ends at k (10)
            model.add_linear_constraint(out_k - in_k == -term[k], name=f"Eq8_10_flow_{kid[k]}")
        # "one and only one arc leaving k" (text of the paper): the route is a path
        model.add_linear_constraint(expr=out_k, ub=1, name=f"OutDeg_{kid[k]}")
        model.add_linear_constraint(expr=in_k, ub=1, name=f"InDeg_{kid[k]}")
    model.add_linear_constraint(mathopt.fast_sum(term.values()) <= p_max, name="Eq11_max_terminals")

    # ---- arrival times along the route
    for k, l in arcs:
        t = t_min[(k, l)]
        if k == gate:
            model.add_linear_constraint(T[l] >= t * x[(k, l)], name=f"Tgate_{kid[l]}")
        else:
            # x_kl = 1: T_l >= T_k + t_kl.   x_kl = 0: vacuous, because T lies in [0, t_cap]
            model.add_linear_constraint(T[k] - T[l] + (t_cap + t) * x[(k, l)] <= t_cap,
                                        name=f"T_{kid[k]}_{kid[l]}")

    # ---- lower level: every destination chooses one transfer stop on the route (AoN)
    for j in J:
        model.add_linear_constraint(mathopt.fast_sum(z[(k, jj)] for k, jj in pairs if jj == j) == 1,
                                    name=f"AoN_{jid[j]}")
    for k in on_route:
        model.add_linear_constraint(on_route[k] == mathopt.fast_sum(x[a] for a in in_arcs[k]),
                                    name=f"OnRouteDef_{kid[k]}")
    for k, j in pairs:
        if k != gate:                                  # users can always get off at the gate
            model.add_linear_constraint(z[(k, j)] <= on_route[k], name=f"OnRoute_{kid[k]}_{jid[j]}")
            # w_kj >= T_k when z_kj = 1
            model.add_linear_constraint(w[(k, j)] >= T[k] - t_cap + t_cap * z[(k, j)],
                                        name=f"Wlink_{kid[k]}_{jid[j]}")

    # ---- route length [NOT in the paper]: every used arc reaches one new stop
    if n_stops is not None:
        used_arcs = mathopt.fast_sum(x.values())
        model.add_linear_constraint(used_arcs == n_stops if stops_mode == "exact"
                                    else used_arcs <= n_stops, name=f"NStops_{stops_mode}")

    return {"model": model, "X": x, "z": z, "b": buses, "T": T, "K": K, "J": J, "gate": gate,
            "demand": demand, "t_min": t_min, "dist": dist, "t_cap": t_cap,
            "params": dict(freq=freq, t_outside_h=t_outside_h, t_rest_h=t_rest_h, c_d=c_d,
                           c_t=c_t, t0=t0, vot=vot, rho_o=rho_o, rho_u=rho_u)}


def extract_solution(built, result):
    if not result.has_primal_feasible_solution():
        raise RuntimeError(f"No feasible solution found: {result.termination.reason}; "
                           f"{result.termination.detail}")
    values = result.variable_values()
    gate = built["gate"]
    next_stop = {a[0]: a[1] for a, v in built["X"].items() if values[v] > 0.5}
    route = [gate]
    while route[-1] in next_stop and len(route) <= len(built["K"]):
        route.append(next_stop[route[-1]])

    # the reported costs are recomputed from the route (independent of the solver's T, w, z, b)
    ev = evaluate_route(route, gate, built["J"], built["demand"], built["t_min"], built["dist"],
                        **built["params"])
    solver_obj = result.objective_value()
    bound = result.best_objective_bound()
    gap = ((solver_obj - bound) / abs(solver_obj)
           if solver_obj and abs(bound) != float("inf") else None)
    return {
        "status": result.termination.reason.name,
        "proven_optimal": result.termination.reason == mathopt.TerminationReason.OPTIMAL,
        "route": route,
        "arrival_min": [round(ev["cum"][k], 2) for k in route],
        "transfer_stop": ev["transfer_stop"],
        "buses": ev["buses"],
        "time_in_city_h": ev["time_in_city_h"],
        "length_km": ev["length_km"],
        "operator_cost": ev["operator_cost"],
        "user_cost": ev["user_cost"],
        "objective": ev["objective"],
        "solver_objective": solver_obj,
        "best_bound": bound,
        "gap": gap,
    }


def load_data(root: Path):
    node_files = [root / "Nodes/N-Bus-Sol.csv", root / "Nodes/N-FGC-Sol.csv",
                  root / "Nodes/N-Metro-Sol.csv", root / "Nodes/N-Tram-Sol.csv"]
    nodes = pd.concat([pd.read_csv(path) for path in node_files], ignore_index=True)
    stops = nodes["id"].tolist()
    if len(stops) != len(set(stops)):
        raise ValueError("Duplicate candidate stop IDs found in node files")

    gate = "Glories - City Gate"
    K = [gate] + stops
    costs = pd.read_csv(root / "Data/Feed_optimization/Cost_matrix.csv")
    demand_table = pd.read_csv(root / "Data/Feed_optimization/Demand_matrix.csv", index_col=0)
    if not {"origen", "dest", "cost", "distance"} <= set(costs.columns):
        raise ValueError("Cost_matrix.csv must contain origen, dest, cost, distance")
    t_min = {(r.origen, r.dest): r.cost for r in costs.itertuples()}
    dist = {(r.origen, r.dest): r.distance for r in costs.itertuples() if pd.notna(r.distance)}
    if demand_table.empty:
        raise ValueError("Demand_matrix.csv contains no demand rows")
    demand_row = demand_table.iloc[0]
    destinations = sorted(set(demand_table.columns) & set(costs.loc[costs["origen"].isin(K), "dest"]))
    if not destinations:
        raise ValueError("No demand destinations overlap the cost matrix")
    demand = {j: float(demand_row[j]) for j in destinations}
    return K, destinations, demand, t_min, dist, gate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).parent)
    parser.add_argument("--solver", choices=("highs", "scip"), default="highs")
    parser.add_argument("--time-limit", type=float, default=1200)
    parser.add_argument("--gap", type=float, default=0.01)
    parser.add_argument("--n-stops", type=int, default=2,
                        help="stops after the gate; omit for free length (paper)")
    parser.add_argument("--stops-mode", choices=("max", "exact"), default="exact")
    parser.add_argument("--t-max", type=float, default=60,
                        help="cap on the minutes inside the city, one way (NOT in the paper)")
    parser.add_argument("--demand-scale", type=float, default=1.0)
    args = parser.parse_args()



    K, J, demand, t_min, dist, gate = load_data(args.root)


    demand = {j: v * args.demand_scale for j, v in demand.items()}
    built = build_model(K, J, gate, demand, t_min, dist,
                        freq=4.4, t_outside_h=45 / 60, t_rest_h=5.8333333 / 60,
                        c_d=0.55, c_t=55.81, t0=5, vot=12, p_max=1,
                        n_stops=args.n_stops, stops_mode=args.stops_mode, t_max_min=args.t_max)
    solver_type = {"highs": mathopt.SolverType.HIGHS, "scip": mathopt.SolverType.GSCIP}[args.solver]
    params = mathopt.SolveParameters(time_limit=timedelta(seconds=args.time_limit),
                                     relative_gap_tolerance=args.gap, enable_output=True)
    result = mathopt.solve(built["model"], solver_type, params=params)

    solution = extract_solution(built, result)
    solution.update(n_stops=args.n_stops, stops_mode=args.stops_mode, t_max_min=built["t_cap"],
                    demand_scale=args.demand_scale, solver=args.solver,
                    solve_time_seconds=result.solve_time().total_seconds())
    for key, value in solution.items():
        print(f"{key}: {value}")

    out_dir = args.root / "Opt_Res"
    out_dir.mkdir(exist_ok=True)
    name = f"solution_cum_{args.n_stops}_stops_x{args.demand_scale:g}.json"
    with open(out_dir / name, "w") as f:
        json.dump(solution, f, indent=2)


if __name__ == "__main__":
    main()