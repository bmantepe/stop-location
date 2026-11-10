"""
Exact enumeration of every route with exactly n_stops stops (1, 2 or 3) after the gate, for the
single-route BNDP of Roca-Riu et al. (2012) with the in-bus time accumulated along the route
(the same cost model as interurban_bnd_ortools.py).

For n_stops <= 3 this gives the guaranteed optimum without any MILP:
    routes = (#gate arcs) * (arcs per stop)^(n_stops-1)      e.g. ~1e3, ~2.5e5, ~6e7 for n = 1, 2, 3
Every route is evaluated with numpy: buses by Eq. (3), cost by Eqs. (1), (6), (7), and each
destination takes the best stop of the route (all-or-nothing).
"""

from __future__ import annotations

import math

import numpy as np

from a5_interurban_bnd_ortools import evaluate_route

BIG = 1e9


def enumerate_best(K, J, gate, demand, t_min, dist, *, n_stops, freq, t_outside_h, t_rest_h,
                   c_d, c_t, t0, vot, rho_o=0.5, rho_u=0.5, t_max_min=None):
    """Returns (best_objective, best_route). t_max_min caps the minutes inside the city."""
    if n_stops not in (1, 2, 3):
        raise ValueError("n_stops must be 1, 2 or 3")
    stops = [k for k in K if k != gate]
    sid = {k: i for i, k in enumerate(stops)}
    N = len(stops)
    dem = np.array([demand[j] for j in J], dtype=float)

    t_g = np.full(N, np.inf)            # gate -> stop arc time
    d_g = np.zeros(N)
    for k in stops:
        if (gate, k) in t_min and (gate, k) in dist:
            t_g[sid[k]] = t_min[(gate, k)]
            d_g[sid[k]] = dist[(gate, k)]
    tau = np.full((N, len(J)), BIG)     # local chain time stop -> destination
    for k in stops:
        for jx, j in enumerate(J):
            if (k, j) in t_min:
                tau[sid[k], jx] = t_min[(k, j)]
    tau_g = np.array([t_min.get((gate, j), BIG) for j in J])

    nb_idx, nb_t, nb_d = [], [], []
    for k in stops:
        cand = [l for l in stops if l != k and (k, l) in t_min and (k, l) in dist]
        nb_idx.append(np.array([sid[l] for l in cand], dtype=int))
        nb_t.append(np.array([t_min[(k, l)] for l in cand], dtype=float))
        nb_d.append(np.array([dist[(k, l)] for l in cand], dtype=float))

    cap = np.inf if t_max_min is None else t_max_min
    fixed_h = t_outside_h + t_rest_h

    def objective(elapsed_min, length_km, user_hours):
        buses = np.ceil((fixed_h + elapsed_min / 60.0) * freq - 1e-9)
        return rho_o * (buses * c_t * t0 + c_d * length_km) + rho_u * vot * user_hours

    best = [math.inf, None]

    def consider(obj, mask, route_of):
        obj = np.where(mask, obj, np.inf)
        i = int(np.argmin(obj))
        if obj[i] < best[0] - 1e-12:
            best[0], best[1] = float(obj[i]), route_of(i)

    for a in range(N):
        if not np.isfinite(t_g[a]) or t_g[a] > cap:
            continue
        cum_a = t_g[a]
        m_a = np.minimum(tau_g, cum_a + tau[a])                 # best stop so far for each destination
        if n_stops == 1:
            obj = objective(cum_a, d_g[a], float(m_a @ dem) / 60.0)
            consider(np.array([obj]), np.array([True]), lambda i, a=a: [gate, stops[a]])
            continue
        b_idx = nb_idx[a]
        if len(b_idx) == 0:
            continue
        cum_b = cum_a + nb_t[a]
        len_b = d_g[a] + nb_d[a]
        if n_stops == 2:
            m_b = np.minimum(m_a[None, :], cum_b[:, None] + tau[b_idx])
            obj = objective(cum_b, len_b, (m_b @ dem) / 60.0)
            consider(obj, cum_b <= cap, lambda i, a=a, b_idx=b_idx: [gate, stops[a], stops[b_idx[i]]])
            continue
        for ib, b in enumerate(b_idx):                          # n_stops == 3
            if cum_b[ib] > cap:
                continue
            c_idx = nb_idx[b]
            if len(c_idx) == 0:
                continue
            cum_c = cum_b[ib] + nb_t[b]
            len_c = len_b[ib] + nb_d[b]
            m_b = np.minimum(m_a, cum_b[ib] + tau[b])
            m_c = np.minimum(m_b[None, :], cum_c[:, None] + tau[c_idx])
            obj = objective(cum_c, len_c, (m_c @ dem) / 60.0)
            consider(obj, (cum_c <= cap) & (c_idx != a),
                     lambda i, a=a, b=b, c_idx=c_idx: [gate, stops[a], stops[b], stops[c_idx[i]]])

    if best[1] is None:
        return math.inf, None
    # recompute the winner with the reference evaluation used by the MILP code
    ev = evaluate_route(best[1], gate, J, demand, t_min, dist, freq=freq, t_outside_h=t_outside_h,
                        t_rest_h=t_rest_h, c_d=c_d, c_t=c_t, t0=t0, vot=vot, rho_o=rho_o, rho_u=rho_u)
    return ev["objective"], best[1]