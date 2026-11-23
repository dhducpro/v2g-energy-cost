"""Individual LPs for the sign of guaranteed site energy-cost saving.

exposures(model, prices) returns dimensionless scores in the native EV order.
prices is one positive purchase price per interval. Time-varying power
availability is supported through the original model's connection matrix.
The classification assumes no aggregate site limit and positive background.
"""
import numpy as np
from scipy.optimize import linprog


def _score(model, vehicle, prices):
    mask = np.asarray(model["availability"][vehicle], dtype=bool)
    duration = np.asarray(model["duration"], dtype=float)[mask]
    prices = prices[mask]
    slots = duration.size
    charge_limit = duration * model["Pc"][vehicle]
    discharge_limit = duration * model["Pd"][vehicle]
    initial = model["E0"][vehicle]
    lower = model["Emin"][vehicle] - initial
    upper = model["Emax"][vehicle] - initial
    target = model["Etarget"][vehicle] - initial
    eta_c = model["eta_c"][vehicle]
    eta_d = model["eta_d"][vehicle]

    if not lower <= 0 <= upper or not 0 < eta_c <= 1 or not 0 < eta_d <= 1:
        raise ValueError("Invalid energy bound or efficiency for " + str(vehicle))
    required = max(target / eta_c, 0.0)
    if required > charge_limit.sum() + 1e-8 or target > upper + 1e-8:
        raise ValueError("Vehicle cannot reach its target by charging alone.")
    if slots == 0 or np.all(discharge_limit == 0):
        return 0.0
    if target == upper:
        # Retaining all charging after withdrawal leaves no room above target.
        return 0.0

    # Under this fleet model the charging-only target is a lower bound, the
    # import price is positive and there are no intermediate mandatory targets.
    # Filling the least expensive connected slots computes C_i^0 exactly.
    baseline = 0.0
    for slot in np.argsort(prices):
        charged = min(required, charge_limit[slot])
        baseline += charged * prices[slot]
        required -= charged
        if required <= 0:
            break

    # Variables are grid-side charged and discharged energy in each connected
    # interval, scaled by the reciprocal of priced discharge; final index is
    # the nonnegative Charnes-Cooper scale.
    eye = np.eye(slots)
    zeros = np.zeros((slots, slots))
    cumulative = np.tril(np.ones((slots, slots)))
    soc = np.column_stack((eta_c * cumulative, -cumulative / eta_d))
    rows = np.vstack((
        np.column_stack((eye, zeros, -charge_limit)),
        np.column_stack((zeros, eye, -discharge_limit)),
        np.column_stack((soc, np.full(slots, -upper))),
        np.column_stack((-soc, np.full(slots, lower))),
        np.r_[-soc[-1], target],
        np.r_[np.full(slots, eta_c), np.zeros(slots), -upper],
    ))
    objective = np.r_[prices, np.zeros(slots), -baseline]
    normalization = np.r_[np.zeros(slots), prices, 0.0]
    result = linprog(
        objective, A_ub=rows, b_ub=np.zeros(len(rows)),
        A_eq=normalization[None, :], b_eq=[1.0], bounds=(0.0, None),
        method="highs", options={"primal_feasibility_tolerance": 1e-9,
                                 "dual_feasibility_tolerance": 1e-9},
    )
    if result.status == 2:
        return 0.0
    if not result.success:
        raise RuntimeError(f"Vehicle {vehicle}: {result.message}")
    if result.fun < -1e-7:
        raise ArithmeticError("Negative incremental charging cost per dollar of discharge.")
    return float(np.clip(1.0 - result.fun, 0.0, 1.0))


def exposures(model, prices):
    """Return the per-EV best positive fraction of gross discharge value retained."""
    prices = np.asarray(prices, dtype=float)
    if prices.shape != (model["T"],) or not np.all(np.isfinite(prices)) or np.any(prices <= 0):
        raise ValueError("Require a finite positive purchase price in every interval.")
    return [_score(model, vehicle, prices) for vehicle in range(model["N"])]
