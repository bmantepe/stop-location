"""Run the MathOpt solver and heuristic-level experiments."""

from __future__ import annotations

import json
from datetime import timedelta
from pathlib import Path

from ortools.math_opt.python import mathopt

from a1_position_mathopt import build_model, extract_solution, load_data


def solve_for_n_stops(
    *,
    n_stops: int,
    data: tuple,
    solver_type: mathopt.SolverType,
    heuristic: mathopt.Emphasis | None,
) -> dict:
    K, J, demand, t_min, dist, gate = data
    built = build_model(
        K,
        J,
        gate,
        demand,
        t_min,
        dist,
        freq=4.4,
        t_outside_h=45 / 60,
        t_rest_h=5.8333333 / 60,
        c_d=0.55,
        c_t=55.81,
        t0=5,
        vot=12,
        p_max=1,
        n_stops=n_stops,
    )
    result = mathopt.solve(
        built["model"],
        solver_type,
        params=mathopt.SolveParameters(
            time_limit=timedelta(hours=1),
            relative_gap_tolerance=0.01,
            heuristics=heuristic,
            enable_output=True,
        ),
    )
    solution = extract_solution(built, result)
    solution["n_stops"] = n_stops
    solution["solve_time_seconds"] = result.solve_time().total_seconds()
    return solution


def main() -> None:
    root = Path(__file__).parent
    output_dir = root / "Opt_Res"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "exact_solutions.json"
    data = load_data(root)

    heuristic_levels = [
        None,
        # mathopt.Emphasis.OFF,
        # mathopt.Emphasis.LOW,
        # mathopt.Emphasis.MEDIUM,
        # mathopt.Emphasis.HIGH,
        # mathopt.Emphasis.VERY_HIGH,
    ]
    solvers = {
        "highs": mathopt.SolverType.HIGHS,
#        "scip": mathopt.SolverType.GSCIP,
    }
    if output_path.exists():
        with output_path.open() as file:
            saved_data = json.load(file)
        results = saved_data.get("results", [])
    else:
        results = []

    for solver_name, solver_type in solvers.items():
        for heuristic in heuristic_levels:
            heuristic_name = heuristic.name if heuristic is not None else None
            for n_stops in range(2, 11):
                # if heuristic is None and solver_name == "highs" and n_stops < 8:
                #     continue
                print(
                    f"Solving solver={solver_name}, "
                    f"heuristics={heuristic_name}, n_stops={n_stops} ...",
                    flush=True,
                )
                solution = solve_for_n_stops(
                    n_stops=n_stops,
                    data=data,
                    solver_type=solver_type,
                    heuristic=heuristic,
                )
                solution["solver"] = solver_name
                solution["heuristics"] = heuristic_name
                solution["n_stops"] = n_stops
                results.append(solution)
                # individual_path = (
                #     output_dir
                #     / f"solution_{solver_name}_{heuristic_name or 'default'}_"
                #     f"{n_stops}_stops.json"
                # )
                # with individual_path.open("w") as file:
                #     json.dump(solution, file, indent=2)
                print(
                    f"  status={solution['status']} "
                    f"objective={solution['objective']} "
                    f"time={solution['solve_time_seconds']:.3f}s",
                    flush=True,
                )
                with output_path.open("w") as file:
                    json.dump({"results": results}, file, indent=2)

if __name__ == "__main__":
    main()
