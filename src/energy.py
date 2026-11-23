"""Scenario LPs and independent schedule audits for fixed charging under discharge withdrawals."""
from datetime import date
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack


def _model(model):
    """Validate and copy inputs; never mutate a supplied case."""
    a = np.asarray(model["availability"], dtype=float).copy()
    if a.ndim != 2 or not a.size:
        raise ValueError("availability must be a nonempty N by T matrix")
    n, nt = a.shape
    if np.any((a != 0) & (a != 1)):
        raise ValueError("Availability must be binary on every interval")
    if int(model.get("N", n)) != n or int(model.get("T", nt)) != nt:
        raise ValueError("N/T do not match availability")
    m = {"name": str(model.get("name", "unnamed")), "N": n, "T": nt,
         "availability": a}
    for key in ("duration", "L", "Ulo", "Uhi"):
        if key not in model:
            raise ValueError(f"Missing explicit model field {key}")
        m[key] = np.broadcast_to(np.asarray(model[key], float), (nt,)).copy()
    for key in ("E0", "Emin", "Emax", "Etarget", "Pc", "Pd", "eta_c", "eta_d"):
        m[key] = np.broadcast_to(np.asarray(model[key], float), (n,)).copy()
    if any(not np.isfinite(v).all() for v in m.values() if isinstance(v, np.ndarray)):
        raise ValueError("Model arrays must be finite")
    if np.any(m["duration"] <= 0):
        raise ValueError("Intervals must have strictly positive duration")
    if np.any(m["Pc"] < 0) or np.any(m["Pd"] < 0):
        raise ValueError("Power ratings must be nonnegative")
    for key in ("eta_c", "eta_d"):
        if np.any(m[key] <= 0) or np.any(m[key] > 1):
            raise ValueError("Efficiencies must be in (0,1]")
    if np.any(m["Ulo"] < 0) or np.any(m["Uhi"] < m["Ulo"]):
        raise ValueError("Coefficient endpoints must obey 0 <= Ulo <= Uhi")
    if (np.any(m["Emax"] < m["Emin"]) or np.any(m["E0"] < m["Emin"] - 1e-10)
            or np.any(m["E0"] > m["Emax"] + 1e-10)
            or np.any(m["Etarget"] < m["Emin"] - 1e-10)
            or np.any(m["Etarget"] > m["Emax"] + 1e-10)):
        raise ValueError("Initial/target energy must lie within battery bounds")
    return m


def _q(q, n):
    if isinstance(q, bool) or int(q) != q or not 0 <= q <= n:
        raise ValueError("q must be an integer from zero through N")
    return int(q)


def audit_schedule(model, q, charge, discharge, M, guard=True):
    """Check physical rows directly, independently of the sparse LP matrix."""
    m = _model(model)
    q = _q(q, m["N"])
    c, d = np.asarray(charge, float), np.asarray(discharge, float)
    if c.shape != m["availability"].shape or d.shape != c.shape:
        raise ValueError("Schedule shape differs from availability")
    if not np.isfinite(c).all() or not np.isfinite(d).all() or not np.isfinite(float(M)):
        raise ValueError("Schedule powers and capacity must be finite")
    dt = m["duration"][None, :]
    energy = m["E0"][:, None] + np.cumsum(
        dt * (m["eta_c"][:, None] * c - d / m["eta_d"][:, None]), axis=1)
    nd = m["E0"][:, None] + np.cumsum(dt * m["eta_c"][:, None] * c, axis=1)
    energy_errors = [0., float(np.max(m["Emin"][:, None] - energy)),
                     float(np.max(energy - m["Emax"][:, None])),
                     float(np.max(m["Etarget"] - energy[:, -1]))]
    if guard:
        energy_errors.append(float(np.max(nd - m["Emax"][:, None])))
    removed = np.sort(d, axis=0)[-q:].sum(axis=0) if q else np.zeros(m["T"])
    net = c.sum(axis=0) - d.sum(axis=0) + removed
    rows = np.maximum(m["L"] + m["Ulo"] * net, m["L"] + m["Uhi"] * net)
    power_error = max(0., -float(M), float(np.max(-c)), float(np.max(-d)),
                      float(np.max(c - m["availability"] * m["Pc"][:, None])),
                      float(np.max(d - m["availability"] * m["Pd"][:, None])),
                      float(np.max(rows - M)))
    overlap = max(0., float(np.max(np.minimum(c, d))))
    return {"energy_residual_kwh": max(energy_errors),
            "power_residual_kw": power_error, "overlap_kw": overlap,
            "guard_checked": bool(guard)}


def prices_for_day(model):
    day = date.fromisoformat(model["metadata"]["horizon_date"])
    edges = np.asarray(model["boundaries_seconds"])
    hours = (edges[:-1] + edges[1:]) / 7200
    off, mid, peak = (.0562, .0925, .267) if 6 <= day.month <= 9 else (.061, .075, .087)
    if day.weekday() >= 5:
        return np.full(hours.size, off)
    return np.select([(hours >= 12) & (hours < 18), (hours >= 8) & (hours < 23)], [peak, mid], off)


class Rows:
    def __init__(self):
        self.values, self.cols, self.idx, self.rhs = [], [], [], []

    def add(self, terms, rhs):
        idx = len(self.rhs)
        self.rhs.append(float(rhs))
        for col, val in terms:
            if val:
                self.idx.append(idx)
                self.cols.append(col)
                self.values.append(val)

    def matrix(self, ncol):
        return coo_matrix((self.values, (self.idx, self.cols)), shape=(len(self.rhs), ncol)).tocsr()


def optimize(model, mode, price_vector=None, site_cap=None, audit=audit_schedule):
    charging_only = mode == "charge_only"
    q = 0 if charging_only else 1
    a = np.asarray(model["availability"])
    n, nt = a.shape
    duration, load = np.asarray(model["duration"]), np.asarray(model["L"])
    ec, ed = np.asarray(model["eta_c"]), np.asarray(model["eta_d"])
    p = prices_for_day(model) if price_vector is None else np.asarray(price_vector)
    cost = duration * p
    row, col = np.nonzero(a)
    ix = {pair: j for j, pair in enumerate(zip(row, col))}
    size = len(row)
    c = lambda i, t: ix[i, t]
    d = lambda i, t: size + ix[i, t]
    g = lambda t: 2 * size + t
    gw = lambda i, t: 2 * size + nt + ix[i, t]
    z = 2 * size + nt + (0 if charging_only else size)
    M = z+1
    nv = M+1
    ub = Rows()
    bounds = [(0, model["Pc"][i]) for i in row]
    bounds.extend((0, 0 if charging_only else model["Pd"][i]) for i in row)
    bounds.extend([(0, None)] * (nv - len(bounds)))
    bounds[M] = (0, site_cap)
    for i in range(n):
        energy, guard = [], []
        for t in np.flatnonzero(a[i]):
            guard.append((c(i,t), duration[t] * ec[i]))
            energy.extend([(c(i,t), duration[t] * ec[i]), (d(i,t), -duration[t] / ed[i])])
            ub.add(energy, model["Emax"][i] - model["E0"][i])
            ub.add([(j, -v) for j,v in energy], model["E0"][i] - model["Emin"][i])
        ub.add([(j, -v) for j,v in energy], model["E0"][i] - model["Etarget"][i])
        ub.add(guard, model["Emax"][i] - model["E0"][i])
    for t in range(nt):
        active = np.flatnonzero(a[:,t])
        net = [(c(i,t), 1) for i in active] + [(d(i,t), -1) for i in active]
        possible = [None] if charging_only else [None, *active]
        for i in possible:
            demand = net if i is None else net + [(d(i,t), 1)]
            meter = g(t) if i is None else gw(i,t)
            ub.add(demand + [(meter, -1)], -load[t])
            for u in {model["Ulo"][t], model["Uhi"][t]}:
                ub.add([(j, u*v) for j,v in demand] + [(M,-1)], -load[t])
    nominal_terms = [(g(t), cost[t]) for t in range(nt)]
    obj = np.zeros(nv)
    for j, v in nominal_terms:
        obj[j] = v
    ub.add(nominal_terms + [(z,-1)], 0)
    if not charging_only:
        for i in range(n):
            scenario = [(gw(i,t) if a[i,t] else g(t), cost[t]) for t in range(nt)]
            ub.add(scenario + [(z,-1)], 0)
    A, b = ub.matrix(nv), np.asarray(ub.rhs)
    options = dict(dual_feasibility_tolerance=1e-9, primal_feasibility_tolerance=1e-9)
    def solve(objective, matrix, rhs):
        result = linprog(objective, A_ub=matrix, b_ub=rhs, bounds=bounds, method="highs", options=options)
        if not result.success:
            raise RuntimeError((model["name"], result.message))
        return result
    edge, cap = np.zeros(nv), np.zeros(nv)
    edge[z], cap[M] = 1, 1
    objectives = {"nominal":obj, "worst":edge, "capacity":cap}
    order = {"charge_only":["nominal", "capacity"],
             "nominal_first":["nominal", "worst", "capacity"],
             "robust_first":["worst", "nominal", "capacity"]}[mode]
    optima = {}
    for name in order:
        objective = objectives[name]
        last = solve(objective, A, b)
        optima[name] = last.fun
        if name != "capacity":
            A = vstack((A, objective.reshape(1,-1)),format="csr")
            b = np.append(b, last.fun + 1e-8)
    if charging_only:
        optima["worst"] = optima["nominal"]
    charge, discharge = np.zeros(a.shape), np.zeros(a.shape)
    charge[row,col] = last.x[:size]
    discharge[row,col] = last.x[size:2*size]
    # Removing a simultaneous charge/discharge cycle preserves the battery
    # trajectory on delivery and leaves more room if the driver withdraws.
    round_trip_efficiency = (ec * ed)[:, None]
    cycle = np.minimum(charge, discharge / round_trip_efficiency)
    charge -= cycle
    discharge -= round_trip_efficiency * cycle

    nominal_import = load + charge.sum(axis=0) - discharge.sum(axis=0)
    nominal_bill = np.maximum(nominal_import, 0) @ cost
    # Each row cancels one driver's discharge for the entire session.
    withdrawal_bills = np.maximum(nominal_import[None, :] + discharge, 0) @ cost
    bill_increases = withdrawal_bills - nominal_bill
    largest_increase = bill_increases.max()
    worst_bill = nominal_bill + max(0, largest_increase)
    driver = int(bill_increases.argmax()) if largest_increase > 1e-7 else None

    peak_ev_load = nominal_import - load + np.max(discharge, axis=0)
    rating = np.maximum(load + model["Ulo"] * peak_ev_load,
                        load + model["Uhi"] * peak_ev_load).max()
    checks = audit(model, q, charge, discharge, optima["capacity"])
    discrepancy = max(abs(nominal_bill - optima["nominal"]), abs(worst_bill - optima["worst"]))
    residuals = ("energy_residual_kwh", "power_residual_kw", "overlap_kw")
    assert discrepancy < 1e-7, (model["name"], discrepancy)
    assert abs(rating - optima["capacity"]) < 1e-6, (model["name"], rating, optima["capacity"])
    assert max(checks[key] for key in residuals) < 1e-6, (model["name"], checks)

    return {
        "nominal_usd": float(nominal_bill),
        "worst_usd": float(worst_bill),
        "max_increment_usd": float(worst_bill - nominal_bill),
        "protected_peak_kw": float(rating),
        "min_net_kw": float(nominal_import.min()),
        "discharge_kwh": float(np.sum(discharge * duration)),
        "withdrawal_ev_1based": driver + 1 if driver is not None else "",
        "withdrawn_discharge_kwh": float(discharge[driver] @ duration) if driver is not None else 0.0,
        "energy_residual_kwh": checks["energy_residual_kwh"],
        "power_residual_kw": checks["power_residual_kw"],
        "overlap_kw": checks["overlap_kw"],
        "economic_check_usd": float(discrepancy),
    }
