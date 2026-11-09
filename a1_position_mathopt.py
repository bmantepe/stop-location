"""Position-based exact formulation of the interurban bus network problem."""

from __future__ import annotations

import argparse
import json
from datetime import timedelta
from pathlib import Path
import pandas as pd

from ortools.math_opt.python import mathopt


def load_data(root: Path):
    node_files = [
        root / "Nodes/N-Bus-Sol.csv",
        root / "Nodes/N-FGC-Sol.csv",
        root / "Nodes/N-Metro-Sol.csv",
        root / "Nodes/N-Tram-Sol.csv",
    ]
    nodes = pd.concat([pd.read_csv(path) for path in node_files], ignore_index=True)
    stops = nodes["id"].tolist()
    if len(stops) != len(set(stops)):
        raise ValueError("Duplicate candidate stop IDs found in node files")

    gate = "Glories - City Gate"
    K = [gate] + stops
    costs = pd.read_csv(root / "Data/Feed_optimization/Cost_matrix.csv")
    demand_table = pd.read_csv(root / "Data/Feed_optimization/Demand_matrix.csv",index_col=0)

    if not {"origen", "dest", "cost", "distance"} <= set(costs.columns):
        raise ValueError("Cost_matrix.csv must contain origen, dest, cost, distance")
    t_min = {(row.origen, row.dest): row.cost for row in costs.itertuples()}
    dist = {
        (row.origen, row.dest): row.distance
        for row in costs.itertuples()
        if pd.notna(row.distance)
    }

    if demand_table.empty:
        raise ValueError("Demand_matrix.csv contains no demand rows")

    demand_row = demand_table.iloc[0]
    destinations = sorted(
        set(demand_table.columns)
        & set(costs.loc[costs["origen"].isin(K), "dest"])
    )
    if not destinations:
        raise ValueError("No demand destinations overlap the cost matrix")
    demand = {j: float(demand_row[j]) for j in destinations}
    return K, destinations, demand, t_min, dist, gate




def build_model(
    K,
    J,
    gate,
    demand,
    t_min,
    dist,
    *,
    n_stops,
    freq,
    t_outside_h,
    t_rest_h,
    c_d,
    c_t,
    t0,
    vot,
    rho_o=0.5,
    rho_u=0.5,
    p_max = 1
):
    """Build an exact-stop position formulation."""
    if n_stops < 1:
        raise ValueError("n_stops must be at least one")

    candidates = [k for k in K if k != gate]
    arcs = [
        (k, l)
        for k in candidates
        for l in candidates
        if k != l and (k, l) in t_min and (k, l) in dist
    ]
    first_arcs = [
        k for k in candidates if (gate, k) in t_min and (gate, k) in dist
    ]
    pairs = [(k, j) for k in K for j in J if (k, j) in t_min]
    uncovered = set(J) - {j for _, j in pairs}
    if uncovered:
        raise ValueError(
            f"Destinations without any reachable transfer stop: {sorted(uncovered)}"
        )

    model = mathopt.Model(name="Interurban_Bus_Network_Position")
    positions = range(n_stops)

    # y[p,k] selects the candidate stop at route position p.
    y = {
        (p, k): model.add_binary_variable(name=f"Y_{p}_{k}")
        for p in positions
        for k in candidates
    }

    # w[p,k,l] selects the transition from position p to p + 1.
    w = {
        (p, k, l): model.add_binary_variable(name=f"W_{p}_{k}_{l}")
        for p in range(n_stops - 1)
        for k, l in arcs
    }

    # z[k,j] assigns destination j to transfer stop k.
    jid = {j: i for i, j in enumerate(J)}
    kid = {k: i for i, k in enumerate(K)}
    z = {
        (k, j): model.add_binary_variable(name=f"Z_{kid[k]}_{jid[j]}")
        for k, j in pairs
    }

    for p in positions:
        model.add_linear_constraint(
            sum(y[p, k] for k in candidates) == 1,
            name=f"OneStopAtPosition_{p}",
        )
    for k in candidates:
        model.add_linear_constraint(
            sum(y[p, k] for p in positions) <= 1,
            name=f"StopUsedOnce_{kid[k]}",
        )
        if k not in first_arcs:
            model.add_linear_constraint(
                y[0, k] == 0,
                name=f"FirstArcUnavailable_{kid[k]}",
            )

    for p in range(n_stops - 1):
        model.add_linear_constraint(
            sum(w[p, k, l] for k, l in arcs)
            == 1,
            name=f"OneTransitionAtPosition_{p}",
        )
        for k in candidates:
            model.add_linear_constraint(
                sum(w[p, k, l] for kk, l in arcs if kk == k)
                == y[p, k],
                name=f"TransitionLeaves_{p}_{kid[k]}",
            )
        for l in candidates:
            model.add_linear_constraint(
                sum(w[p, k, ll] for k, ll in arcs if ll == l)
                == y[p + 1, l],
                name=f"TransitionEnters_{p}_{kid[l]}",
            )

    # Cumulative arrival time at each route position, in hours.
    max_arc_time_h = max(
        (t_min[a] / 60.0 for a in [(gate, k) for k in first_arcs] + arcs),
        default=0.0,
    )
    max_route_time_h = n_stops * max_arc_time_h
    arrival = {
        p: model.add_variable(
            lb=0.0,
            ub=max_route_time_h,
            name=f"ArrivalPosition_{p}",
        )
        for p in positions
    }
    model.add_linear_constraint(
        arrival[0]
        == sum((t_min[gate, k] / 60.0) * y[0, k] for k in first_arcs),
        name="ArrivalAtFirstStop",
    )
    for p in range(n_stops - 1):
        model.add_linear_constraint(
            arrival[p + 1]
            == arrival[p]
            + sum((t_min[k, l] / 60.0) * w[p, k, l] for k, l in arcs),
            name=f"ArrivalRecurrence_{p}",
        )

    # Arrival time indexed by stop, zero when a stop is not selected.
    arrival_stop = {
        k: model.add_variable(
            lb=0.0,
            ub=max_route_time_h,
            name=f"ArrivalStop_{kid[k]}",
        )
        for k in candidates
    }
    for k in candidates:
        for p in positions:
            model.add_linear_constraint(
                arrival_stop[k] <= arrival[p] + max_route_time_h * (1 - y[p, k]),
                name=f"ArrivalStopUpper_{p}_{kid[k]}",
            )
            model.add_linear_constraint(
                arrival_stop[k] >= arrival[p] - max_route_time_h * (1 - y[p, k]),
                name=f"ArrivalStopLower_{p}_{kid[k]}",
            )
        model.add_linear_constraint(
            arrival_stop[k]
            <= max_route_time_h * sum(y[p, k] for p in positions),
            name=f"ArrivalStopUnused_{kid[k]}",
        )

    # Transfer-time products arrival_stop[k] * z[k,j].
    transfer_time = {
        (k, j): model.add_variable(
            lb=0.0,
            ub=max_route_time_h,
            name=f"TransferTime_{kid[k]}_{jid[j]}",
        )
        for k, j in pairs
        if k != gate
    }
    for k, j in pairs:
        if k == gate:
            continue
        model.add_linear_constraint(
            transfer_time[k, j] <= arrival_stop[k],
            name=f"TransferTimeArrivalUpper_{kid[k]}_{jid[j]}",
        )
        model.add_linear_constraint(
            transfer_time[k, j] <= max_route_time_h * z[k, j],
            name=f"TransferTimeZUpper_{kid[k]}_{jid[j]}",
        )
        model.add_linear_constraint(
            transfer_time[k, j]
            >= arrival_stop[k] - max_route_time_h * (1 - z[k, j]),
            name=f"TransferTimeLower_{kid[k]}_{jid[j]}",
        )

    for j in J:
        model.add_linear_constraint(
            sum(z[k, jj] for k, jj in pairs if jj == j) == 1,
            name=f"AoN_{jid[j]}",
        )
    for k in candidates:
        on_route = sum(y[p, k] for p in positions)
        for kk, j in pairs:
            if kk == k:
                model.add_linear_constraint(
                    z[k, j] <= on_route,
                    name=f"OnRoute_{kid[k]}_{jid[j]}",
                )

    route_time_h = arrival[n_stops - 1]
    route_distance_km = sum(
        dist[gate, k] / 1.0 * y[0, k] for k in first_arcs
    ) + sum(dist[k, l] * w[p, k, l] for p, k, l in w)
    buses = model.add_integer_variable(lb=0, name="Buses")
    operator_cost = buses * c_t * t0 + c_d * route_distance_km
    model.add_linear_constraint(
        buses >= (t_outside_h + route_time_h + t_rest_h) * freq,
        name="Eq3_buses",
    )

    user_cost = vot * sum(
        demand[j]
        * (
            (
                transfer_time[k, j]
                if k != gate
                else 0.0
            )
            + (t_min[k, j] / 60.0) * z[k, j]
        )
        for k, j in pairs
    )
    model.minimize(rho_o * operator_cost + rho_u * user_cost)

    return {
        "model": model,
        "y": y,
        "w": w,
        "z": z,
        "b": buses,
        "arrival": arrival,
        "arrival_stop": arrival_stop,
        "transfer_time": transfer_time,
        "route_time_h": route_time_h,
        "route_distance_km": route_distance_km,
        "operator_cost": operator_cost,
        "user_cost": user_cost,
        "K": K,
        "J": J,
        "candidates": candidates,
        "gate": gate,
        "kid": kid,
        "jid": jid,
        "demand": demand,
        "t_min": t_min,
        "dist": dist,
        "vot": vot,
        "c_d": c_d,
        "c_t": c_t,
        "t0": t0,
        "rho_o": rho_o,
        "rho_u": rho_u,
    }


def extract_solution(built, result):
    if not result.has_primal_feasible_solution():
        raise RuntimeError(
            f"No feasible solution found: {result.termination.reason}; "
            f"{result.termination.detail}"
        )
    values = result.variable_values()
    route = [built["gate"]]
    for p in range(len(built["arrival"])):
        selected = [
            k for k in built["candidates"] if values[built["y"][p, k]] > 0.5
        ]
        route.extend(selected)

    transfer_stop = {
        j: k
        for (k, j), variable in built["z"].items()
        if values[variable] > 0.5
    }
    route_time_h = values[built["route_time_h"]]
    distance_km = sum(
        built["dist"][built["gate"], k] * values[built["y"][0, k]]
        for k in built["candidates"]
        if (built["gate"], k) in built["dist"]
    ) + sum(
        built["dist"][k, l] * values[variable]
        for (p, k, l), variable in built["w"].items()
    )
    operator_cost = (
        values[built["b"]] * built["c_t"] * built["t0"]
        + built["c_d"] * distance_km
    )
    user_cost = built["vot"] * sum(
        built["demand"][j]
        * (
            (
                values[built["transfer_time"][k, j]]
                if k != built["gate"]
                else 0.0
            )
            + built["t_min"][k, j] / 60.0 * values[variable]
        )
        for (k, j), variable in built["z"].items()
    )
    objective = built["rho_o"] * operator_cost + built["rho_u"] * user_cost
    return {
        "status": result.termination.reason.name,
        "proven_optimal": result.termination.reason == mathopt.TerminationReason.OPTIMAL,
        "route": route,
        "transfer_stop": transfer_stop,
        "buses": round(values[built["b"]]),
        "time_in_city_h": route_time_h,
        "length_km": distance_km,
        "operator_cost": operator_cost,
        "user_cost": user_cost,
        "objective": objective,
        "best_bound": result.best_objective_bound(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).parent)
    parser.add_argument("--n-stops", type=int, default=2)
    parser.add_argument("--time-limit", type=float, default=60)
    parser.add_argument("--gap", type=float, default=0.01)
    args = parser.parse_args()

    K, J, demand, t_min, dist, gate = load_data(args.root)
    built = build_model(
        K,
        J,
        gate,
        demand,
        t_min,
        dist,
        n_stops=args.n_stops,
        freq=4.4,
        t_outside_h=45 / 60,
        t_rest_h=5.8333333 / 60,
        c_d=0.55,
        c_t=55.81,
        t0=5,
        vot=12,
    )
    result = mathopt.solve(
        built["model"],
        mathopt.SolverType.HIGHS,
        params=mathopt.SolveParameters(
            time_limit=timedelta(seconds=args.time_limit),
            relative_gap_tolerance=args.gap,
            enable_output=True,
        ),
    )
    solution = extract_solution(built, result)
    solution["n_stops"] = args.n_stops
    solution["solve_time_seconds"] = result.solve_time().total_seconds()
    for key, value in solution.items():
        print(f"{key}: {value}")

    output_dir = args.root / "Opt_Res"
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / f"solution_position_{args.n_stops}_stops.json").open("w") as file:
        json.dump(solution, file, indent=2)


if __name__ == "__main__":
    main()
