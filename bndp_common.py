"""Shared data loading, parameters and route evaluation for the interurban BNDP.

Implements the cost functions of Roca-Riu, Estrada & Trapote (2012),
"The design of interurban bus networks in city centers", for a single
interurban route (Mataró -> Barcelona, line e11) entering Barcelona through
the city gate.

Units: times in minutes in the data, converted to hours for costs; distances in km.
"""

from __future__ import annotations

import heapq
import math
import random
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
COST_MATRIX = ROOT / "Data/Feed_optimization/Cost_matrix.csv"
DEMAND_MATRIX = ROOT / "Data/Feed_optimization/Demand_matrix.csv"
CANDIDATE_FILES = [ROOT / f"Nodes/N-{m}-Sol.csv" for m in ("Bus", "FGC", "Tram", "Metro")]
POI_FILE = ROOT / "Nodes/N-POIs.csv"
GATE_ID = "Glories - City Gate"


@dataclass
class Params:
    rho_operator: float = 0.5
    rho_user: float = 0.5
    vot: float = 12.0                # €/h
    p_terminals: int = 1             # max number of terminals |P| <= p
    frequency: float = 4.4           # buses/h
    t_outside_h: float = 45 / 60     # T_A, origin -> gate (h)
    t_rest_h: float = 5.8333 / 60    # T_rest (h)
    cd: float = 0.55                 # €/km
    ct: float = 55.81                # €/(bus·h)
    t0: float = 5.0                  # hours of service
    w_iu: float = 1.0                # perception weight of in-city interurban in-vehicle time


@dataclass
class Instance:
    stops: list[str]                                   # candidate stops (excluding gate)
    pois: list[str]
    demand: dict[str, float]                           # poi -> trips
    arcs: dict[tuple[str, str], tuple[float, float]]   # (from, to) -> (minutes, km); from may be the gate
    poi_cost: dict[str, dict[str, float]]              # stop -> {poi: perceived minutes stop -> poi}
    names: dict[str, str] = field(default_factory=dict)
    poi_names: dict[str, str] = field(default_factory=dict)

    def out_arcs(self) -> dict[str, list[tuple[str, float, float]]]:
        out: dict[str, list[tuple[str, float, float]]] = {}
        for (k, l), (t, d) in self.arcs.items():
            out.setdefault(k, []).append((l, t, d))
        return out


def load_instance(
    demand_row: str | None = None, sample_stops: int | None = None, seed: int = 0
) -> Instance:
    """Load the full instance, or a random subset of `sample_stops` candidate stops (for testing)."""
    cm = pd.read_csv(COST_MATRIX, dtype={"origen": str, "dest": str}).drop_duplicates()
    dups = cm.duplicated(["origen", "dest"]).sum()
    if dups:
        raise ValueError(f"{dups} (origen, dest) pairs with conflicting values in the cost matrix")

    poi_rows = cm[cm["distance"].isna()]
    arc_rows = cm[cm["distance"].notna()]

    dm = pd.read_csv(DEMAND_MATRIX, index_col=0)
    row = dm.loc[demand_row] if demand_row else dm.iloc[0]
    demand = {str(c): float(v) for c, v in row.items() if v > 0}

    cand = pd.read_csv('Nodes/N-All-Sol.csv', dtype={"id": str})
    names = dict(zip(cand["id"], cand["name"]))
    names[GATE_ID] = GATE_ID
    poi_df = pd.read_csv(POI_FILE, dtype={"id": str})
    poi_names = dict(zip(poi_df["id"], poi_df["id"]))

    # A stop is usable only if it is a candidate, has POI costs and can be reached by the bus.
    stops_with_poi = set(poi_rows["origen"])
    arc_nodes = set(arc_rows["dest"])
    stops = sorted(set(cand["id"]) & stops_with_poi & arc_nodes)
    stop_set = set(stops)

    arcs = {}
    for k, l, t, d in arc_rows[["origen", "dest", "cost", "distance"]].itertuples(index=False):
        if (k == GATE_ID or k in stop_set) and l in stop_set and k != l:
            arcs[(k, l)] = (float(t), float(d))

    poi_cost: dict[str, dict[str, float]] = {k: {} for k in stops}
    for k, j, c in poi_rows[["origen", "dest", "cost"]].itertuples(index=False):
        if k in stop_set and j in demand:
            poi_cost[k][j] = float(c)

    pois = sorted(demand, key=str)
    missing = [j for j in pois if not any(j in poi_cost[k] for k in stops)]
    if missing:
        raise ValueError(f"POIs without any stop -> POI path: {missing}")

    inst = Instance(stops=stops, pois=pois, demand=demand, arcs=arcs, poi_cost=poi_cost, names=names,
                    poi_names=poi_names)
    return subset_instance(inst, sample_stops, seed) if sample_stops else inst


def subset_instance(inst: Instance, n_stops: int, seed: int = 0) -> Instance:
    """Random subset of `n_stops` candidate stops (arcs and POI costs restricted accordingly)."""
    if n_stops >= len(inst.stops):
        return inst
    stops = sorted(random.Random(seed).sample(inst.stops, n_stops))
    s = set(stops)
    arcs = {(k, l): v for (k, l), v in inst.arcs.items() if l in s and (k == GATE_ID or k in s)}
    return Instance(
        stops=stops,
        pois=inst.pois,
        demand=inst.demand,
        arcs=arcs,
        poi_cost={k: inst.poi_cost[k] for k in stops},
        names=inst.names,
        poi_names=inst.poi_names,
    )


def hops_from_gate(inst: Instance) -> dict[str, int]:
    out = inst.out_arcs()
    hops = {GATE_ID: 0}
    q = deque([GATE_ID])
    while q:
        k = q.popleft()
        for l, _, _ in out.get(k, []):
            if l not in hops:
                hops[l] = hops[k] + 1
                q.append(l)
    return hops


def shortest_time_from_gate(inst: Instance) -> dict[str, float]:
    out = inst.out_arcs()
    dist = {GATE_ID: 0.0}
    pq = [(0.0, GATE_ID)]
    while pq:
        dk, k = heapq.heappop(pq)
        if dk > dist.get(k, math.inf):
            continue
        for l, t, _ in out.get(k, []):
            nd = dk + t
            if nd < dist.get(l, math.inf):
                dist[l] = nd
                heapq.heappush(pq, (nd, l))
    return dist


def evaluate_route(inst: Instance, prm: Params, route: list[str]) -> dict:
    """Exact cost of a route (sequence of stops after the gate), independent of any solver.

    Lower level (All-or-Nothing): each POI demand alights at the stop of the route that
    minimises in-vehicle time from the gate + perceived time stop -> POI.
    """
    if not route:
        raise ValueError("empty route")
    if len(set(route)) != len(route):
        raise ValueError("route visits a stop twice")

    arrival, t_city, length = {}, 0.0, 0.0
    prev = GATE_ID
    for k in route:
        if (prev, k) not in inst.arcs:
            raise ValueError(f"arc {prev} -> {k} does not exist")
        t, d = inst.arcs[(prev, k)]
        t_city += t
        length += d
        arrival[k] = t_city
        prev = k

    assignment, user_minutes, drops = {}, 0.0, []
    for j in inst.pois:
        best = min(
            ((prm.w_iu * arrival[k] + inst.poi_cost[k][j], k) for k in route if j in inst.poi_cost[k]),
            default=None,
        )
        if best is None:
            raise ValueError(f"POI {j} unreachable from route")
        k = best[1]
        assignment[j] = k
        user_minutes += inst.demand[j] * best[0]
        drops.append({
            "poi_id": j,
            "poi_name": inst.poi_names.get(j, j),
            "demand": inst.demand[j],
            "drop_stop_id": k,
            "drop_stop_name": inst.names.get(k, k),
            "in_vehicle_min": arrival[k],
            "stop_to_poi_min": inst.poi_cost[k][j],
            "total_perceived_min": best[0],
        })

    buses = math.ceil((prm.t_outside_h + t_city / 60 + prm.t_rest_h) * prm.frequency - 1e-9)
    c_op = buses * prm.ct * prm.t0 + prm.cd * length
    c_user = prm.vot * user_minutes / 60
    return {
        "route": list(route),
        "n_stops": len(route),
        "city_time_min": t_city,
        "length_km": length,
        "buses": buses,
        "operator_cost": c_op,
        "user_cost": c_user,
        "objective": prm.rho_operator * c_op + prm.rho_user * c_user,
        "avg_user_time_min": user_minutes / sum(inst.demand.values()),
        "arrival_min": arrival,
        "assignment": assignment,
        "poi_drop_stops": drops,
    }
