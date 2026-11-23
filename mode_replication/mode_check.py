"""Check the uncapped one-withdrawal energy model with explicit binary EV modes."""
import argparse
import csv
import json
import platform
import time
import warnings
from datetime import date
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix


class Variables:
    def __init__(self):
        self.lower = []
        self.upper = []
        self.integer = []

    def add(self, n, lower=0, upper=np.inf, integer=False):
        start = len(self.lower)
        self.lower.extend(np.broadcast_to(lower, n))
        self.upper.extend(np.broadcast_to(upper, n))
        self.integer.extend([int(integer)] * n)
        return np.arange(start, start+n)


class Constraints:
    def __init__(self):
        self.rows = []
        self.cols = []
        self.values = []
        self.lower = []
        self.upper = []

    def add(self, cols, values, low=-np.inf, high=np.inf):
        self.rows.extend([len(self.lower)] * len(cols))
        self.cols.extend(cols)
        self.values.extend(values)
        self.lower.append(low)
        self.upper.append(high)

    def compile(self, n):
        matrix = coo_matrix((self.values, (self.rows, self.cols)), shape=(len(self.lower), n)).tocsc()
        return LinearConstraint(matrix, self.lower, self.upper)


def tariff(model):
    day = date.fromisoformat(model["metadata"]["horizon_date"])
    boundaries = np.asarray(model["boundaries_seconds"])
    h = (boundaries[1:] + boundaries[:-1]) / 7200
    off, mid, high = (.0562, .0925, .267) if 6 <= day.month <= 9 else (.061, .075, .087)
    if day.weekday() >= 5:
        return np.full(len(h), off)
    return np.select([(12 <= h) & (h < 18), (8 <= h) & (h < 23)], [high, mid], default=off)


def control(model, binary, limit, prices=None):
    available = np.array(model["availability"], bool)
    n, nt = available.shape
    i, t = np.nonzero(available)
    nactive = len(i)
    dt = np.asarray(model["duration"])
    price_weight = (tariff(model) if prices is None else np.asarray(prices, float)) * dt
    load = np.asarray(model["L"])
    pc, pd = np.asarray(model["Pc"]), np.asarray(model["Pd"])
    ec, ed = np.asarray(model["eta_c"]), np.asarray(model["eta_d"])
    e0, emin, emax, target = (np.asarray(model[x]) for x in ("E0", "Emin", "Emax", "Etarget"))
    v, con = Variables(), Constraints()
    charge = v.add(nactive, upper=pc[i])
    discharge = v.add(nactive, upper=pd[i])
    state = v.add(nactive, lower=emin[i], upper=emax[i])
    mode = v.add(nactive, upper=1, integer=True) if binary else None
    net = v.add(nt, lower=-np.inf)
    meter = v.add((n+1) * nt).reshape(n+1, nt)
    z = v.add(1)[0]
    active = np.full((n, nt), -1, int)
    active[i, t] = np.arange(nactive)
    for ev in range(n):
        ids = np.flatnonzero(i == ev)
        for index, cell in enumerate(ids):
            cols = [state[cell], charge[cell], discharge[cell]]
            coeff = [1, -dt[t[cell]]*ec[ev], dt[t[cell]]/ed[ev]]
            rhs = e0[ev]
            if index:
                cols.append(state[ids[index-1]])
                coeff.append(-1)
                rhs = 0
            con.add(cols, coeff, low=rhs, high=rhs)
            if binary:
                con.add([charge[cell], mode[cell]], [1, -pc[ev]], high=0)
                con.add([discharge[cell], mode[cell]], [1, pd[ev]], high=pd[ev])
        con.add([state[ids[-1]]], [1], low=target[ev])
        con.add(charge[ids], dt[t[ids]]*ec[ev], high=emax[ev]-e0[ev])
    for slot in range(nt):
        ids = np.flatnonzero(t == slot)
        con.add(np.r_[net[slot], charge[ids], discharge[ids]],
                np.r_[1, -np.ones(len(ids)), np.ones(len(ids))], low=load[slot], high=load[slot])
        for scenario in range(n+1):
            cols, coefs = [meter[scenario, slot], net[slot]], [1, -1]
            if scenario and available[scenario-1, slot]:
                cols.append(discharge[active[scenario-1, slot]])
                coefs.append(-1)
            con.add(cols, coefs, low=0)
    for g in meter:
        con.add(np.r_[g, z], np.r_[price_weight, -1], high=0)
    objective = np.zeros(len(v.lower))
    objective[z] = 1
    constraints = con.compile(len(objective))
    bounds = Bounds(v.lower, v.upper)
    # SciPy forwards these HiGHS controls unchanged.
    options = dict(time_limit=limit, mip_rel_gap=0.0, threads=1,
                   mip_feasibility_tolerance=1e-8, primal_feasibility_tolerance=1e-9)
    start = time.perf_counter()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unrecognized options detected.*")
        sol = milp(objective, integrality=v.integer if binary else None, bounds=bounds,
                   constraints=constraints, options=options)
    elapsed = time.perf_counter()-start
    result = dict(status=sol.status, solver_seconds=elapsed, objective_usd=sol.fun,
                  mip_gap=sol.get("mip_gap"), binary_count=nactive if binary else 0)
    if sol.x is None:
        return result
    c, d = np.zeros((n, nt)), np.zeros((n, nt))
    c[i,t], d[i,t] = sol.x[charge], sol.x[discharge]
    nominal = load+np.sum(c-d, axis=0)
    meters = np.maximum(np.vstack([nominal, nominal+d]), 0)
    result["delivered_min_net_kw"] = np.min(nominal)
    result["delivered_export_kwh"] = np.maximum(-nominal, 0)@dt
    result["worst_bill_usd"] = (meters@price_weight).max()
    energy = e0[:,None]+np.cumsum((c*ec[:,None]-d/ed[:,None])*dt, axis=1)
    canceled_energy = e0[:,None]+np.cumsum(c*ec[:,None]*dt, axis=1)
    result["battery_violation_kwh"] = max(0, np.max(emin[:,None]-energy),np.max(energy-emax[:,None]),
                np.max(canceled_energy-emax[:,None]), np.max(target-energy[:,-1]))
    result["overlap_kw"] = max(0, np.minimum(c,d).max())
    if sol.success:
        assert abs(result["objective_usd"]-result["worst_bill_usd"]) < 1e-6
        assert result["battery_violation_kwh"] < 1e-6
        if binary:
            assert result["overlap_kw"] < 1e-6
    return result


def main():
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--fleets", type=Path, required=True)
    arg.add_argument("--published", type=Path, required=True, help="Existing q_sweep.csv")
    arg.add_argument("--out", type=Path, default=Path("mode_checks.csv"))
    arg.add_argument("--pattern", default="*.json")
    arg.add_argument("--limit", type=float, default=20)
    args = arg.parse_args()
    with args.published.open(newline="") as fh:
        reference = {r["case"]:float(r["robust_upper_usd"]) for r in csv.DictReader(fh) if r["q"]=="1"}
    fields = ["case","date","N","T","binary_count","LP_status","MILP_status","LP_usd","MILP_usd",
              "reference_usd","MILP_minus_LP_usd","LP_seconds","MILP_seconds","MIP_gap",
              "MILP_audit_bill_usd","MILP_overlap_kw","MILP_battery_violation_kwh",
              "MILP_delivered_min_net_kw","MILP_delivered_export_kwh"]
    print(platform.platform(), "Python",platform.python_version(),"SciPy",scipy.__version__, flush=True)
    paths = sorted(args.fleets.glob(args.pattern))
    if not paths:
        raise ValueError("No input files matched")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as out:
        writer=csv.DictWriter(out,fieldnames=fields)
        writer.writeheader()
        for path in paths:
            model=json.loads(path.read_text(encoding="utf8"))
            lp, mip = (control(model, b, args.limit) for b in (False, True))
            ref=reference[model["name"]]
            gap = None if mip["objective_usd"] is None else mip["objective_usd"]-lp["objective_usd"]
            assert lp["status"]==0 and abs(lp["objective_usd"]-ref)<1e-6
            if mip["status"]==0:
                assert abs(gap)<1e-6
            writer.writerow(dict(case=model["name"],date=model["metadata"]["horizon_date"],
                N=model["N"],T=model["T"],binary_count=mip["binary_count"],LP_status=lp["status"],
                MILP_status=mip["status"],LP_usd=lp["objective_usd"],MILP_usd=mip["objective_usd"],
                reference_usd=ref,MILP_minus_LP_usd=gap,LP_seconds=lp["solver_seconds"],
                MILP_seconds=mip["solver_seconds"],MIP_gap=mip["mip_gap"],
                MILP_audit_bill_usd=mip.get("worst_bill_usd"),MILP_overlap_kw=mip.get("overlap_kw"),
                MILP_battery_violation_kwh=mip.get("battery_violation_kwh"),
                MILP_delivered_min_net_kw=mip.get("delivered_min_net_kw"),
                MILP_delivered_export_kwh=mip.get("delivered_export_kwh")))
            out.flush()
            print(path.stem,mip["status"],round(mip["solver_seconds"],3),gap,flush=True)


if __name__=="__main__":
    main()
