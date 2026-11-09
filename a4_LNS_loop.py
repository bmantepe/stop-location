"""Solve the CP-SAT LNS model for routes with 2 through 10 stops."""

from __future__ import annotations

import json
from datetime import timedelta
from pathlib import Path

from ortools.sat import sat_parameters_pb2
from ortools.math_opt.python import mathopt

from a3_LNS_position_cpsat import build_model, extract_solution, load_data


def solve_for_n_stops(data, n_stops: int, *, time_limit: float = 3600) -> dict:
    """Solve one exact-stop-count instance with CP-SAT LNS."""
    K, J, demand, t_sec, dist_m, gate = data
    built = build_model(
        K,
        J,
        gate,
        demand,
        t_sec,
        dist_m,
        freq=4.4,
        t_outside_sec=2700,
        t_rest_sec=350,
        c_d=0.55,
        c_t=55.81,
        t0=5,
        vot=12,
        demand_scale=1000,
        p_max=1,
        n_stops=n_stops,
        stops_mode="exact",
    )

    params = mathopt.SolveParameters(
        time_limit=timedelta(seconds=time_limit),
        cp_sat=sat_parameters_pb2.SatParameters(
            num_search_workers=8,
            log_search_progress=True,
        ),
        relative_gap_tolerance=0.01,
        enable_output=True,
    )
    result = mathopt.solve(
        built["model"],
        mathopt.SolverType.CP_SAT,
        params=params,
    )

    solution = extract_solution(built, result)
    solution.update(
        {
            "solver": "CP_SAT",
            "n_stops": n_stops,
            "stops_mode": "exact",
            "time_limit_seconds": time_limit,
            "solve_time_seconds": result.solve_time().total_seconds(),
            "num_search_workers": 8,
            "relative_gap": 0.01,
        }
    )
    return solution


def main() -> None:
    root = Path(__file__).parent
    output_dir = root / "Opt_Res"
    output_dir.mkdir(parents=True, exist_ok=True)

    data = load_data(root, demand_scale=1000)
    results = []

    for n_stops in range(2, 11):
        print(f"Solving CP-SAT LNS for n_stops={n_stops} ...", flush=True)
        solution = solve_for_n_stops(data, n_stops)
        results.append(solution)

        individual_path = output_dir / f"solution_cpsat_{n_stops}_stops.json"
        with individual_path.open("w") as file:
            json.dump(solution, file, indent=2)

        print(
            f"  status={solution['status']} "
            f"objective={solution['objective']} "
            f"time={solution['solve_time_seconds']:.3f}s",
            flush=True,
        )

        with (output_dir / "solutions_cpsat_2_to_10_stops.json").open("w") as file:
            json.dump({"results": results}, file, indent=2)


if __name__ == "__main__":
    main()
