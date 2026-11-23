"""Direct clipped-meter cost optimization on the temporal Caltech fleets.

Usage: python validate_threshold.py --fleets PATH/TO/observed_cases --out PATH
Requires Python 3.10+, NumPy, SciPy (HiGHS) and the adjacent criterion.py.
Whole-session withdrawals, zero credit for negative metered import, one common
nonadaptive charge schedule, positive energy rates and no exogenous site cap.
"""
import argparse
import csv
import json
from datetime import date
from pathlib import Path
from itertools import combinations
from math import comb
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp
from scipy.sparse import coo_matrix, csr_matrix, hstack, vstack


def prices_for_day(model):
    day = date.fromisoformat(model["metadata"]["horizon_date"])
    e = np.asarray(model["boundaries_seconds"])
    h = (e[:-1] + e[1:]) / 7200
    off, mid, peak = (.0562, .0925, .267) if 6 <= day.month <= 9 else (.061, .075, .087)
    if day.weekday() >= 5:
        return np.full(h.size, off)
    return np.select([(h >= 12) & (h < 18), (h >= 8) & (h < 23)], [peak, mid], off)


class Rows:
    def __init__(self):
        self.r, self.c, self.v, self.b = [], [], [], []

    def add(self, cols, values, rhs):
        self.r.extend([len(self.b)] * len(cols))
        self.c.extend(cols)
        self.v.extend(values)
        self.b.append(rhs)

    def build(self, ncol):
        return coo_matrix((self.v, (self.r, self.c)), shape=(len(self.b), ncol)).tocsr(), np.array(self.b)


class Fleet:
    def __init__(self, model):
        self.model = model
        self.n, self.t = model["N"], model["T"]
        self.availability = np.asarray(model["availability"])
        i, t = np.nonzero(self.availability)
        self.i, self.slot = i, t
        self.nz = len(i)
        self.dt = np.array(model["duration"])
        self.load = np.array(model["L"])
        self.p = prices_for_day(model)
        self.w = self.dt * self.p
        self.ec, self.ed = map(np.array, (model["eta_c"], model["eta_d"]))
        self.e0, self.emin, self.emax, self.target = map(np.array, (model["E0"], model["Emin"], model["Emax"], model["Etarget"]))
        self.z, self.theta = 2*self.nz, 2*self.nz+1
        self.r = self.theta + 1 + np.arange(self.n)
        self.ncols = int(self.r[-1]+1)
        self.bounds = ([(0,j) for j in np.array(model["Pc"])[i]] + [(0,j) for j in np.array(model["Pd"])[i]]
                       + [(0,None)]*(self.ncols-2*self.nz))
        self.bounds = list(self.bounds)
        rows = Rows()
        for j in range(self.n):
            cols = np.flatnonzero(i==j)
            ct, ds = self.dt[t[cols]] * self.ec[j], -self.dt[t[cols]] / self.ed[j]
            for h in range(1,len(cols)+1):
                sel = np.r_[cols[:h], self.nz + cols[:h]]
                vals = np.r_[ct[:h],ds[:h]]
                rows.add(sel, vals, self.emax[j]-self.e0[j])
                lower = self.target[j] if h==len(cols) else self.emin[j]
                rows.add(sel,-vals,self.e0[j]-max(lower, self.emin[j]))
            rows.add(cols,ct,self.emax[j]-self.e0[j])
            rows.add(np.r_[self.nz+cols,self.theta,self.r[j]],np.r_[self.w[t[cols]],-1,-1],0)
        self.A,self.b = rows.build(self.ncols)
        self.obj = np.zeros(self.ncols)
        self.obj[self.z] = 1

    def recover(self,x):
        c,d = np.zeros((self.n,self.t)),np.zeros((self.n,self.t))
        c[self.i,self.slot] = x[:self.nz]
        d[self.i,self.slot] = x[self.nz:2*self.nz]
        k = (self.ec*self.ed)[:,None]
        cycles = np.minimum(c,d/k)
        c -= cycles
        d -= k*cycles
        energy = self.e0[:,None]+np.cumsum((c*self.ec[:,None]-d/self.ed[:,None])*self.dt,axis=1)
        guard = self.e0+np.sum(c*self.ec[:,None]*self.dt,axis=1)
        epserr = max(0, np.max(energy-self.emax[:,None]),np.max(self.emin[:,None]-energy),np.max(self.target-energy[:,-1]),np.max(guard-self.emax))
        pwerr = max(0,np.max(np.minimum(c,d)),np.max(-c),np.max(-d),
                    np.max(c-self.availability*np.array(self.model["Pc"])[:,None]),
                    np.max(d-self.availability*np.array(self.model["Pd"])[:,None]))
        assert epserr < 1e-6 and pwerr < 1e-6
        return self.load + np.sum(c-d,axis=0), d, float(epserr), float(pwerr)

    def solve(self, A, b, bounds=None):
        opt = linprog(self.obj,A_ub=A,b_ub=b,bounds=self.bounds if bounds is None else bounds,method="highs",
                      options={"primal_feasibility_tolerance":1e-9,"dual_feasibility_tolerance":1e-9})
        if not opt.success:
            raise RuntimeError((self.model["name"],opt.message))
        return opt

    def scenario_row(self, owners, active_hours):
        """Bound the purchase bill for one set of withdrawals and import intervals."""
        charging = np.flatnonzero(active_hours[self.slot])
        discharging = charging[~np.isin(self.i[charging], owners)]

        columns = np.r_[charging, self.nz + discharging, self.z]
        weights = np.r_[self.w[self.slot[charging]], -self.w[self.slot[discharging]], -1.]
        row = csr_matrix(
            (weights, (np.zeros(len(weights), dtype=int), columns)),
            shape=(1, self.ncols),
        )
        background_bill = self.w[active_hours] @ self.load[active_hours]
        return row, -float(background_bill)

    def explicit(self,q):
        # Build a separate epigraph for every possible withdrawn set, using a
        # distinct metered import variable in each interval (no signed shortcut).
        # Metered cost is increasing in removed discharge, so exactly q owners
        # suffice even though the contract allows up to q.
        scenarios = list(combinations(range(self.n),q))
        rows = Rows()
        for s,owners in enumerate(scenarios):
            grid = self.ncols+s*self.t+np.arange(self.t)
            for h in range(self.t):
                slots = np.flatnonzero(self.slot==h)
                dslots = slots[~np.isin(self.i[slots],owners)]
                cols = np.r_[slots,self.nz+dslots,grid[h]]
                vals = np.r_[np.ones(slots.size),-np.ones(dslots.size),-1.]
                rows.add(cols,vals,-self.load[h])
            rows.add(np.r_[grid,self.z],np.r_[self.w,-1.],0.)
        extra = len(scenarios)*self.t
        B,b = rows.build(self.ncols+extra)
        A = vstack((hstack((self.A,csr_matrix((self.A.shape[0],extra)))),B))
        ans = linprog(np.r_[self.obj,np.zeros(extra)],A_ub=A,b_ub=np.r_[self.b,b],
                      bounds=self.bounds+[(0,None)]*extra,method="highs",
                      options={"primal_feasibility_tolerance":1e-9,"dual_feasibility_tolerance":1e-9})
        if not ans.success:
            raise RuntimeError(ans.message)
        net,d,energy,power = self.recover(ans.x)
        bills = [float(np.maximum(net+np.sum(d[list(owners)],axis=0),0) @ self.w) for owners in scenarios]
        gap = abs(max(bills)-ans.fun)
        assert gap<1e-6
        return float(ans.fun), gap, energy, power

    def charge_only(self):
        # With no station cap, choose each driver's cheapest available kWh.
        total=float(self.load @ self.w)
        for j in range(self.n):
            need=max(0.,(self.target[j]-self.e0[j])/self.ec[j])
            for t in np.argsort(self.p,kind="stable"):
                amount=min(need,self.dt[t]*self.availability[j,t]*self.model['Pc'][j])
                total+=amount*self.p[t]
                need-=amount
            assert need<1e-8
        return total

    def _withdrawals_milp(self, net, discharge, q, contributors):
        """Choose the worst withdrawals when too many sets exist to enumerate."""
        supply = discharge[contributors]
        count = len(contributors)
        max_import = net + np.sort(supply, axis=0)[-q:].sum(axis=0)
        importing = net >= 0
        # Only these intervals can switch from export to import after a withdrawal.
        mixed_hours = np.flatnonzero((net < 0) & (max_import > 0))

        # Variables: one withdrawal flag per driver, then an import flag and
        # a metered purchase for each interval that can switch sign.
        num_mixed = len(mixed_hours)
        num_flags = count + num_mixed
        num_vars = num_flags + num_mixed
        drivers = np.arange(count)
        objective = np.r_[
            -(supply[:, importing] @ self.w[importing]),
            np.zeros(num_mixed),
            -self.w[mixed_hours],
        ]
        rows = Rows()
        rows.add(drivers, np.ones(count), q)
        for index, hour in enumerate(mixed_hours):
            flag = count + index
            purchase = num_flags + index
            rows.add([flag, purchase], [-max_import[hour], 1], 0)
            # With the flag on, the purchase cannot exceed import after withdrawals.
            rows.add(np.r_[drivers, flag, purchase], np.r_[-supply[:, hour], -net[hour], 1], 0)

        matrix, upper = rows.build(num_vars)
        result = milp(
            objective,
            integrality=np.r_[np.ones(num_flags), np.zeros(num_mixed)],
            bounds=Bounds(np.zeros(num_vars), np.r_[np.ones(num_flags), max_import[mixed_hours]]),
            constraints=LinearConstraint(matrix, np.full(len(upper), -np.inf), upper),
            options={"mip_rel_gap": 0.},
        )
        if not result.success:
            raise RuntimeError(result.message)
        return tuple(contributors[result.x[:count] > .5])

    def worst(self, net, d, q):
        """Find the highest metered bill and the drivers whose withdrawals cause it."""
        discharge_value = d @ self.w
        contributors = np.flatnonzero(discharge_value > 0)
        if q >= len(contributors):
            candidates = [tuple(contributors)]
        elif q == 0:
            candidates = [()]
        elif min(net) >= -1e-10:
            # The station already imports in every interval; only total discharge value matters.
            candidates = [tuple(np.argsort(-discharge_value)[:q])]
        elif comb(len(contributors), q) <= 30000:
            candidates = combinations(contributors, q)
        else:
            candidates = [self._withdrawals_milp(net, d, q, contributors)]

        worst_bill = -np.inf
        worst_owners = import_hours = None
        for owners in candidates:
            grid_import = net + np.sum(d[list(owners)], axis=0)
            bill = float(np.maximum(grid_import, 0) @ self.w)
            if bill > worst_bill:
                worst_bill = bill
                worst_owners = owners
                import_hours = grid_import > 0
        return worst_bill, worst_owners, import_hours

    def protected(self, q):
        """Add withdrawal scenarios until the worst metered bill meets the LP bound."""
        # Start with the signed bill: crediting exports gives a lower bound on costs.
        row, rhs = self.scenario_row((), np.ones(self.t, bool))
        budget_row = csr_matrix(
            (np.r_[float(q), np.ones(self.n)],
             (np.zeros(self.n + 1, int), np.r_[self.theta, self.r])),
            shape=(1, self.ncols),
        )
        matrix = vstack((self.A, row + budget_row), format="csr")
        bounds = np.r_[self.b, rhs]
        cuts = 0
        while True:
            solution = self.solve(matrix, bounds)
            if cuts == 0:
                signed_lower = float(solution.fun)
            net, discharge, energy_error, power_error = self.recover(solution.x)
            bill, owners, import_hours = self.worst(net, discharge, q)
            gap = bill - solution.fun
            if gap < 1e-7:
                return float(solution.fun), float(bill), cuts, energy_error, power_error, signed_lower

            row, rhs = self.scenario_row(owners, import_hours)
            matrix = vstack((matrix, row), format="csr")
            bounds = np.r_[bounds, rhs]
            cuts += 1
            if cuts > 600:
                raise RuntimeError(("too many cuts", self.model["name"], q, gap))


def write_csv(path,rows):
    with path.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fleets",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    from criterion import exposures
    details,fleets,drivers=[],[],[]
    start=time.perf_counter()
    for path in sorted(args.fleets.glob("*.json")):
        m=json.loads(path.read_text())
        f=Fleet(m)
        # Charge-only is exactly q=N: the all-withdrawn scenario pays for all
        # planned charging, and zero discharge is feasible with minimum charging.
        base,base_ub,*_=f.protected(f.n)
        base_exact=f.charge_only()
        assert max(abs(base-base_exact),abs(base_ub-base_exact))<1e-6
        vals=exposures(m,f.p)
        beta=float(sum(vals))
        cutoff=int(np.ceil(beta-1e-8))
        name=m["name"]
        for i,val in enumerate(vals):
            drivers.append(dict(case=name,driver=i+1,beta=float(val)))
        observed=None
        for q in range(f.n+1):
            lb,ub,cuts,ep,pp,signed=f.protected(q)
            direct,directgap= None,None
            if q<=1 or f.n<=6:
                direct,fg,eu,pu=f.explicit(q)
                directgap=max(abs(direct-lb),abs(direct-ub),fg)
                assert directgap<1e-6
                ep,pp=max(ep,eu),max(pp,pu)
            positive=base-ub>1e-6
            negligible=base_ub-lb<=1e-6
            assert positive or negligible
            if not positive and observed is None:
                observed=q
            agree=positive==(q<cutoff)
            details.append(dict(case=name,date=m["metadata"]["horizon_date"],N=f.n,q=q,beta_sum=beta,predicted_first_zero=cutoff,
                                charge_only_usd=base_exact,robust_lower_usd=lb,robust_upper_usd=ub,saving_usd=base_exact-ub,
                                positive=positive,prediction_matches=agree,cuts=cuts,initial_signed_lower_usd=signed,
                                clipping_cost_gap_usd=ub-signed,primal_dual_gap_usd=ub-lb,
                                explicit_clipped_lp_usd=direct,explicit_check_usd=directgap,energy_residual_kwh=ep,power_residual_kw=pp))
        datestr=m["metadata"]["horizon_date"]
        fleets.append(dict(case=name,date=datestr,N=f.n,beta_sum=beta,first_zero_predicted=cutoff,first_zero_observed=observed,
                           charge_only_greedy_check_usd=abs(base-base_exact)))
        print(f"{name} threshold={observed} beta={beta:.6g} elapsed={time.perf_counter()-start:.1f}",flush=True)
    write_csv(args.out/"q_sweep.csv",details)
    write_csv(args.out/"fleet_thresholds.csv",fleets)
    write_csv(args.out/"driver_exposures.csv",drivers)
    print(f"n={len(details)}, matches={sum(r['prediction_matches'] for r in details)}")


if __name__=="__main__":
    main()
