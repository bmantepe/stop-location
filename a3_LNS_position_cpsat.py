"""Position-based CP-SAT/LNS formulation of the bus network problem."""

from __future__ import annotations

import argparse
import json
import math
from datetime import timedelta
from fractions import Fraction
from pathlib import Path

from ortools.math_opt.python import mathopt
from ortools.sat import sat_parameters_pb2

def load_data(root: Path, demand_scale=1):
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
    demand_table = pd.read_csv(
        root / "Data/Feed_optimization/Demand_matrix.csv", index_col=0
    )

    if not {"origen", "dest", "cost", "distance"} <= set(costs.columns):
        raise ValueError("Cost_matrix.csv must contain origen, dest, cost, distance")

    # The original model treated cost as minutes and distance as kilometres.
    # Convert to integer seconds and metres for CP-SAT.
    t_sec = {
        (row.origen, row.dest): int(round(float(row.cost) * 60))
        for row in costs.itertuples()
        if pd.notna(row.cost)
    }
    dist_m = {
        (row.origen, row.dest): int(round(float(row.distance) * 1000))
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

    # Preserve up to three decimal places of demand while making it integral.
    demand = {
        j: int(round(float(demand_row[j]) * demand_scale))
        for j in destinations
    }
    return K, destinations, demand, t_sec, dist_m, gate

def build_model(
    K,
    J,
    gate,
    demand,
    t_sec,
    dist_m,
    *,
    n_stops,
    freq,
    t_outside_sec,
    t_rest_sec,
    c_d,
    c_t,
    t0,
    vot,
    rho_o=0.5,
    rho_u=0.5,
    demand_scale=1,
):
    candidates = [k for k in K if k != gate]
    arcs = [
        (k, l)
        for k in candidates
        for l in candidates
        if k != l and (k, l) in t_sec and (k, l) in dist_m
    ]
    first_stops = [
        k for k in candidates if (gate, k) in t_sec and (gate, k) in dist_m
    ]
    pairs = [(k, j) for k in K for j in J if (k, j) in t_sec]
    uncovered = set(J) - {j for _, j in pairs}
    if uncovered:
        raise ValueError(f"Unreachable destinations: {sorted(uncovered)}")

    kid = {k: i for i, k in enumerate(K)}
    jid = {j: i for i, j in enumerate(J)}
    model = mathopt.Model(name="Interurban_Bus_Network_Position_CP_SAT")

    y = {
        (p, k): model.add_binary_variable(name=f"Y_{p}_{kid[k]}")
        for p in range(n_stops)
        for k in candidates
    }
    w = {
        (p, k, l): model.add_binary_variable(name=f"W_{p}_{kid[k]}_{kid[l]}")
        for p in range(n_stops - 1)
        for k, l in arcs
    }
    z = {
        (k, j): model.add_binary_variable(name=f"Z_{kid[k]}_{jid[j]}")
        for k, j in pairs
    }

    for p in range(n_stops):
        model.add_linear_constraint(
            sum(y[p, k] for k in candidates) == 1,
            name=f"OneStopAtPosition_{p}",
        )
    for k in candidates:
        model.add_linear_constraint(
            sum(y[p, k] for p in range(n_stops)) <= 1,
            name=f"StopUsedOnce_{kid[k]}",
        )
        if k not in first_stops:
            model.add_linear_constraint(
                y[0, k] == 0,
                name=f"FirstArcUnavailable_{kid[k]}",
            )

    for p in range(n_stops - 1):
        for k in candidates:
            model.add_linear_constraint(
                sum(w[p, k, l] for kk, l in arcs if kk == k) == y[p, k],
                name=f"TransitionLeaves_{p}_{kid[k]}",
            )
        for l in candidates:
            model.add_linear_constraint(
                sum(w[p, k, ll] for k, ll in arcs if ll == l) == y[p + 1, l],
                name=f"TransitionEnters_{p}_{kid[l]}",
            )

    max_arc_time = max(
        [t_sec[gate, k] for k in first_stops]
        + [t_sec[a] for a in arcs],
        default=0,
    )
    max_route_time = n_stops * max_arc_time
    arrival = {
        p: model.add_integer_variable(
            lb=0, ub=max_route_time, name=f"ArrivalPosition_{p}"
        )
        for p in range(n_stops)
    }
    model.add_linear_constraint(
        arrival[0] == sum(t_sec[gate, k] * y[0, k] for k in first_stops),
        name="ArrivalAtFirstStop",
    )
    for p in range(n_stops - 1):
        model.add_linear_constraint(
            arrival[p + 1]
            == arrival[p] + sum(t_sec[k, l] * w[p, k, l] for k, l in arcs),
            name=f"ArrivalRecurrence_{p}",
        )

    arrival_stop = {
        k: model.add_integer_variable(
            lb=0, ub=max_route_time, name=f"ArrivalStop_{kid[k]}"
        )
        for k in candidates
    }
    for k in candidates:
        used = sum(y[p, k] for p in range(n_stops))
        model.add_linear_constraint(
            arrival_stop[k] <= max_route_time * used,
            name=f"ArrivalStopUnused_{kid[k]}",
        )
        for p in range(n_stops):
            model.add_linear_constraint(
                arrival_stop[k] <= arrival[p] + max_route_time * (1 - y[p, k]),
                name=f"ArrivalStopUpper_{p}_{kid[k]}",
            )
            model.add_linear_constraint(
                arrival_stop[k] >= arrival[p] - max_route_time * (1 - y[p, k]),
                name=f"ArrivalStopLower_{p}_{kid[k]}",
            )

    transfer_time = {
        (k, j): model.add_integer_variable(
            lb=0, ub=max_route_time, name=f"TransferTime_{kid[k]}_{jid[j]}"
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
            transfer_time[k, j] <= max_route_time * z[k, j],
            name=f"TransferTimeZUpper_{kid[k]}_{jid[j]}",
        )
        model.add_linear_constraint(
            transfer_time[k, j]
            >= arrival_stop[k] - max_route_time * (1 - z[k, j]),
            name=f"TransferTimeLower_{kid[k]}_{jid[j]}",
        )

    for j in J:
        model.add_linear_constraint(
            sum(z[k, jj] for k, jj in pairs if jj == j) == 1,
            name=f"AoN_{jid[j]}",
        )
    for k in candidates:
        used = sum(y[p, k] for p in range(n_stops))
        for kk, j in pairs:
            if kk == k:
                model.add_linear_constraint(
                    z[k, j] <= used,
                    name=f"OnRoute_{kid[k]}_{jid[j]}",
                )

    route_time = arrival[n_stops - 1]
    route_distance = (
        sum(dist_m[gate, k] * y[0, k] for k in first_stops)
        + sum(dist_m[k, l] * w[p, k, l] for p, k, l in w)
    )
    freq_fraction = Fraction(str(freq))
    buses_ub = (
        freq_fraction.numerator
        * (t_outside_sec + max_route_time + t_rest_sec)
        + 3600 * freq_fraction.denominator
        - 1
    ) // (3600 * freq_fraction.denominator)
    buses = model.add_integer_variable(lb=0, ub=max(0, buses_ub), name="Buses")
    model.add_linear_constraint(
        3600 * freq_fraction.denominator * buses
        >= freq_fraction.numerator
        * (t_outside_sec + route_time + t_rest_sec),
        name="Eq3_buses",
    )

    rho_o_f = Fraction(str(rho_o))
    rho_u_f = Fraction(str(rho_u))
    coefficients = [
        rho_o_f * Fraction(str(c_t)) * Fraction(str(t0)),
        rho_o_f * Fraction(str(c_d)) / 1000,
        rho_u_f * Fraction(str(vot)) / (3600 * demand_scale),
    ]
    objective_scale = math.lcm(*(c.denominator for c in coefficients))
    bus_coef = int(coefficients[0] * objective_scale)
    distance_coef = int(coefficients[1] * objective_scale)
    user_coef = int(coefficients[2] * objective_scale)
    user_seconds = sum(
        (
            demand[j] * transfer_time[k, j]
            if k != gate
            else 0
        )
        + demand[j] * t_sec[k, j] * z[k, j]
        for k, j in pairs
    )
    model.minimize(
        bus_coef * buses + distance_coef * route_distance + user_coef * user_seconds
    )

    return {
        "model": model,
        "y": y,
        "z": z,
        "b": buses,
        "arrival": arrival,
        "transfer_time": transfer_time,
        "route_distance": route_distance,
        "route_time": route_time,
        "objective_scale": objective_scale,
        "K": K,
        "J": J,
        "candidates": candidates,
        "gate": gate,
        "kid": kid,
        "demand": demand,
        "t_sec": t_sec,
        "dist_m": dist_m,
        "demand_scale": demand_scale,
        "c_d": c_d,
        "c_t": c_t,
        "t0": t0,
        "vot": vot,
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
        route.append(
            next(
                k
                for k in built["candidates"]
                if values[built["y"][p, k]] > 0.5
            )
        )
    transfer_stop = {
        j: k
        for (k, j), variable in built["z"].items()
        if values[variable] > 0.5
    }
    time_in_city_h = values[built["route_time"]] / 3600
    length_m = sum(
        built["dist_m"][a, b]
        for a, b in zip(route, route[1:])
    )
    length_km = length_m / 1000
    operator_cost = values[built["b"]] * built["c_t"] * built["t0"]
    operator_cost += built["c_d"] * length_km
    user_cost = built["vot"] * sum(
        (
            (
                values[built["transfer_time"][k, j]]
                if k != built["gate"]
                else 0
            )
            + built["t_sec"][k, j] * values[variable]
        )
        * built["demand"][j]
        for (k, j), variable in built["z"].items()
    ) / (built["demand_scale"] * 3600)
    objective = built["rho_o"] * operator_cost + built["rho_u"] * user_cost
    return {
        "status": result.termination.reason.name,
        "proven_optimal": result.termination.reason == mathopt.TerminationReason.OPTIMAL,
        "route": route,
        "transfer_stop": transfer_stop,
        "buses": round(values[built["b"]]),
        "time_in_city_h": time_in_city_h,
        "length_km": length_km,
        "operator_cost": operator_cost,
        "user_cost": user_cost,
        "objective": objective,
        "best_bound": result.best_objective_bound() / built["objective_scale"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).parent)
    parser.add_argument("--n-stops", type=int, default=2)
    parser.add_argument("--time-limit", type=float, default=7200)
    parser.add_argument("--gap", type=float, default=0.01)
    args = parser.parse_args()

    K, J, demand, t_sec, dist_m, gate = load_data(args.root, demand_scale=1000)
    built = build_model(
        K, J, gate, demand, t_sec, dist_m,
        n_stops=args.n_stops,
        freq=4.4,
        t_outside_sec=2700,
        t_rest_sec=350,
        c_d=0.55,
        c_t=55.81,
        t0=5,
        vot=12,
    )
    result = mathopt.solve(
        built["model"],
        mathopt.SolverType.CP_SAT,
        params=mathopt.SolveParameters(
            time_limit=timedelta(seconds=args.time_limit),
            cp_sat=sat_parameters_pb2.SatParameters(
                num_search_workers=8,
                log_search_progress=True,
            ),
            relative_gap_tolerance=args.gap,
            enable_output=True,
        ),
    )
    solution = extract_solution(built, result)
    solution.update({
        "solver": "CP_SAT",
        "formulation": "position",
        "n_stops": args.n_stops,
        "solve_time_seconds": result.solve_time().total_seconds(),
    })
    for key, value in solution.items():
        print(f"{key}: {value}")

    output_dir = args.root / "Opt_Res"
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / f"solution_cpsat_position_{args.n_stops}_stops.json").open("w") as file:
        json.dump(solution, file, indent=2)


if __name__ == "__main__":
    main()
