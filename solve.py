"""
Interurban Bus Network Design Problem (BNDP) -- Roca-Riu, Estrada & Trapote (2012)
Single-route version, implemented with PuLP.

Paper equations are referenced as (n) in the comments.

Why a single-level MILP is exact here
-------------------------------------
Upper level: min  rho_O*C_O + rho_U*C_U   (7)
Lower level: min  C_U                     (12)
For a fixed route, C_O does not depend on the users' choices and rho_U > 0, so the
upper level also wants the lower-level optimum. Solving both levels at once gives
the same solution as the bilevel formulation.
"""
import pandas as pd
import pulp


def build_model(K, J, GATE, d, t_min, dist, *,
                freq, t_outside_h, t_rest_h, c_d, c_t, t0, vot,
                rho_o=0.5, rho_u=0.5, p_max=1, extras=True,
                n_stops=None, stops_mode="max"):
    """
    K       : list of candidate stops (must include GATE, I = {GATE}, I subset of K)
    J       : list of destination names (centroids)
    d       : dict  destination -> demand d_j  (single origin, so i is dropped)
    t_min   : dict  (a, b) -> time in MINUTES for
                      (stop, stop)    in-vehicle interurban time      t_kl,IU
                      (GATE, stop)    in-vehicle interurban time      t_ik,IU
                      (stop, dest)    local chain: transfer + waiting + urban + egress,
                                      already weighted by w_T, w_W, w_TU, w_E
    dist    : dict  (stop, stop) -> km  (l_kl)
    freq    : buses per hour (f^r)
    t_outside_h, t_rest_h : hours (T_A^r, T_rest)
    c_d [EUR/km], c_t [EUR/h], t0 [h], vot [EUR/h]
    p_max   : maximum number of terminals p (11)
    n_stops : (NOT in the paper) number of stops after the gate. None = free length, as in the paper.
    stops_mode : "max" -> at most n_stops stops, "exact" -> exactly n_stops stops
    extras  : add the constraints that the paper's text implies but does not write
              (in/out degree <= 1 and subtour elimination). Set False for the literal model.
    """
    kid = {k: i for i, k in enumerate(K)}      # integer ids -> safe PuLP variable names
    jid = {j: i for i, j in enumerate(J)}

    # ---------- index sets ----------
    ARCS = [(k, l) for k in K for l in K if k != l and (k, l) in t_min and (k, l) in dist]
    t_ik = {k: (0.0 if k == GATE else t_min.get((GATE, k))) for k in K}   # t_ik,IU
    PAIRS = [(k, j) for k in K for j in J if (k, j) in t_min and t_ik[k] is not None]
    uncovered = set(J) - {j for _, j in PAIRS}
    assert not uncovered, f"Destinations without any reachable transfer stop: {uncovered}"
    out_arcs = {k: [a for a in ARCS if a[0] == k] for k in K}
    in_arcs = {k: [a for a in ARCS if a[1] == k] for k in K}

    # ---------- decision variables ----------
    prob = pulp.LpProblem("Interurban_Bus_Network_Design", pulp.LpMinimize)
    # Upper level: X_kl = 1 if arc (k,l) is used by the route (15)
    X = {a: pulp.LpVariable(f"X_{kid[a[0]]}_{kid[a[1]]}", cat="Binary") for a in ARCS}
    # Terminal indicator (k in P)
    term = {k: pulp.LpVariable(f"Term_{kid[k]}", cat="Binary") for k in K if k != GATE}
    # Lower level: z_kj = 1 if users going to j leave the interurban bus at stop k (transfer stop)
    z = {p: pulp.LpVariable(f"Z_{kid[p[0]]}_{jid[p[1]]}", cat="Binary") for p in PAIRS}
    # b = total number of buses (integer replaces the ceiling in (3))
    b = pulp.LpVariable("Buses", lowBound=0, cat="Integer")

    # ---------- operator cost (1), (2), (4), (5) ----------
    T_C = pulp.lpSum((t_min[a] / 60.0) * X[a] for a in ARCS)      # (4) hours inside the city
    l_r = pulp.lpSum(dist[a] * X[a] for a in ARCS)                # (5) km inside the city
    C_O = b * c_t * t0 + c_d * l_r                                # (1)

    # (3) b >= (T_A + T_C + T_rest) * f.  b has a positive cost, so at the optimum
    #     b equals the smallest integer satisfying this, i.e. the ceiling.
    prob += b >= (t_outside_h + T_C + t_rest_h) * freq, "Eq3_buses"

    # ---------- user cost (6) ----------
    C_U = vot * pulp.lpSum(
        d[j] * ((t_ik[k] + t_min[(k, j)]) / 60.0) * z[(k, j)] for (k, j) in PAIRS
    )

    # ---------- objective (7) ----------
    prob += rho_o * C_O + rho_u * C_U, "Eq7_objective"

    # ---------- upper-level constraints ----------
    for k in K:
        out_k = pulp.lpSum(X[a] for a in out_arcs[k])
        in_k = pulp.lpSum(X[a] for a in in_arcs[k])
        if k == GATE:
            prob += out_k - in_k == 1, "Eq9_route_starts_at_gate"                 # (9)
        else:
            # term_k = 0 -> out - in = 0 (8) intermediate / unused stop
            # term_k = 1 -> in - out = 1 (10) route ends at a terminal
            prob += out_k - in_k == -term[k], f"Eq8_10_flow_{kid[k]}"

    prob += pulp.lpSum(term.values()) <= p_max, "Eq11_max_terminals"             # (11)

    # ---------- lower-level constraints ----------
    # every destination is served through exactly one transfer stop (AoN)
    for j in J:
        prob += pulp.lpSum(z[(k, jj)] for (k, jj) in PAIRS if jj == j) == 1, f"AoN_{jid[j]}"

    # (13)+(14): transfer stops are candidate stops (built into the index sets) and users can
    # only leave where the route stops. The paper states this in words only.
    on_route = {k: pulp.LpVariable(f"On_{kid[k]}", 0, 1) for k in K if k != GATE}
    for k in on_route:
        prob += on_route[k] == pulp.lpSum(X[a] for a in in_arcs[k]), f"OnRoute_def_{kid[k]}"
    for (k, j) in PAIRS:
        if k != GATE:
            prob += z[(k, j)] <= on_route[k], f"OnRoute_{kid[k]}_{jid[j]}"

    # ---------- NOT in the paper's equations (needed for a correct MILP) ----------
    if extras:
        N = len(K)
        for k in K:
            prob += pulp.lpSum(X[a] for a in out_arcs[k]) <= 1, f"OutDeg_{kid[k]}"  # "one and only one arc leaving k"
            prob += pulp.lpSum(X[a] for a in in_arcs[k]) <= 1, f"InDeg_{kid[k]}"
        # MTZ subtour elimination: stops the model from adding a disconnected loop whose
        # stops would count as "on the route" for the lower level
        u = {k: pulp.LpVariable(f"u_{kid[k]}", lowBound=1, upBound=max(N - 1, 1))
             for k in K if k != GATE}
        for (k, l) in ARCS:
            if k != GATE and l != GATE:
                prob += u[k] - u[l] + (N - 1) * X[(k, l)] <= N - 2, f"MTZ_{kid[k]}_{kid[l]}"

    # ---------- OPTIONAL: control the route length (NOT in the paper) ----------
    # Each used arc reaches one new stop, so sum(X) = number of stops after the gate.
    # Requires extras=True (otherwise loops would also be counted).
    if n_stops is not None:
        n_used = pulp.lpSum(X.values())
        if stops_mode == "exact":
            prob += n_used == n_stops, "NStops_exact"
        else:
            prob += n_used <= n_stops, "NStops_max"

    return dict(prob=prob, X=X, z=z, b=b, T_C=T_C, l_r=l_r, C_O=C_O, C_U=C_U)


def extract_solution(model, K, GATE):
    prob = model["prob"]
    # Without this check, a run that stops on the time limit without an integer solution
    # would still print the (fractional) LP values as if they were a route.
    if prob.sol_status not in (pulp.LpSolutionOptimal, pulp.LpSolutionIntegerFeasible):
        raise RuntimeError(f"No integer-feasible solution found (status: {pulp.LpStatus[prob.status]}, "
                           f"sol_status: {prob.sol_status}). Prune arcs, raise timeLimit or set gapRel.")
    X, z = model["X"], model["z"]
    nxt = {a[0]: a[1] for a, v in X.items() if v.value() is not None and v.value() > 0.5}
    route = [GATE]
    while route[-1] in nxt and len(route) <= len(K):
        route.append(nxt[route[-1]])
    transfer = {j: k for (k, j), v in z.items() if v.value() is not None and v.value() > 0.5}
    return dict(
        # PuLP calls any CBC run with an integer solution "Optimal", even when it stopped on the time
        # limit, so the label is taken from sol_status instead.
        status=("Optimal" if prob.sol_status == pulp.LpSolutionOptimal
                else "Feasible, NOT proven optimal (time/gap limit) - check the gap in the CBC log"),
        proven_optimal=(prob.sol_status == pulp.LpSolutionOptimal),
        route=route,
        transfer_stop=transfer,
        buses=int(round(model["b"].value())),
        time_in_city_h=pulp.value(model["T_C"]),
        length_km=pulp.value(model["l_r"]),
        operator_cost=pulp.value(model["C_O"]),
        user_cost=pulp.value(model["C_U"]),
        objective=pulp.value(model["prob"].objective),
    )


# =====================================================================================
if __name__ == "__main__":
    # ---------------- DATA (your loading code, tidied) ----------------
    bus_sol = pd.read_csv('Nodes/N-Bus-Sol.csv')
    fgc_sol = pd.read_csv('Nodes/N-FGC-Sol.csv')
    metro_sol = pd.read_csv('Nodes/N-Metro-Sol.csv')
    tram_sol = pd.read_csv('Nodes/N-Tram-Sol.csv')
    pois = pd.read_csv('Nodes/N-POIs.csv')

    GATE = 'Glories - City Gate'
    stops = pd.concat([bus_sol, fgc_sol, metro_sol, tram_sol], ignore_index=True)['id'].tolist()
    K = [GATE] + stops
    J = pois['poi_name'].tolist()
    assert len(set(K)) == len(K), "duplicate stop ids"

    # ASSUMED FORMATS -- adapt the column names to your files:
    #   N-Costs.csv : origin, destination, time_min, distance_km
    #       (gate->stop, stop->stop, stop->destination). For stop->destination, time_min must
    #       already be the weighted chain  wT*t_T + wW*t_W + wTU*t_TU + wE*t_E  (minutes).
    #   N-Demand.csv: one row, one column per destination (same names as poi_name)
    costs = pd.read_csv('Data/Feed_optimization/Cost_matrix.csv')
    t_min = {(r.origen, r.dest): r.cost for r in costs.itertuples()}
    dist = {(r.origen, r.dest): r.distance for r in costs.itertuples()}
    demand = pd.read_csv('Data/Feed_optimization/Demand_matrix.csv')
    d = demand.iloc[0].to_dict()
    # (long format instead?)  d = dict(zip(demand['poi_name'], demand['demand']))

    # ---------------- CONSTANTS ----------------
    FREQ = 4.4                   # buses per hour
    T_A = 45 / 60                # h, time outside the city
    T_REST = 5.8333333 / 60      # h
    C_D = 0.55                   # EUR/km
    C_T = 55.81                  # EUR/h  <-- NOT in your data: placeholder (paper's value)
    T0 = 5                    # h      <-- NOT in your data: placeholder (paper's 2 h rush hour)
    VOT = 12                     # EUR/h
    RHO_O = RHO_U = 0.5
    N_STOPS = 3            # stops after the gate; use None for free length (paper)
    P_MAX = 1                    # p in (11); with one route there is exactly one terminal

    model = build_model(K, J, GATE, d, t_min, dist,
                        freq=FREQ, t_outside_h=T_A, t_rest_h=T_REST,
                        c_d=C_D, c_t=C_T, t0=T0, vot=VOT,
                        rho_o=RHO_O, rho_u=RHO_U, p_max=P_MAX, extras=True,
                        n_stops=N_STOPS, stops_mode="max")
    # CBC (default). For big models HiGHS is usually much faster:  pip install highspy
    solver = pulp.HiGHS(msg=True, timeLimit=600, gapRel=0.01)
    #solver = pulp.PULP_CBC_CMD(msg=True, timeLimit=600, gapRel=0.01)
    model["prob"].solve(solver)
    sol = extract_solution(model, K, GATE)
    for key, val in sol.items():
        print(f"{key}: {val}")