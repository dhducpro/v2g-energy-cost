"""Recalculate the common physical import limits and all matched cost comparisons."""

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix


def minimum_unidirectional_peak(model):
    """Smallest maximum import serving the supplied targets when discharge is zero."""
    availability = np.asarray(model["availability"], dtype=bool)
    n, slots = availability.shape
    vehicle, interval = np.nonzero(availability)
    duration = np.asarray(model["duration"])
    power = np.asarray(model["Pc"])
    load = np.asarray(model["L"])
    efficiency = np.asarray(model["eta_c"])
    energy = np.maximum(0, np.subtract(model["Etarget"], model["E0"]))
    entries = np.arange(len(vehicle))
    peak = len(vehicle)
    # Charging-only battery energy is nondecreasing. Reaching the departure
    # target exactly respects every prefix bound when arrival and target
    # energies are within their battery limits, as they are in these inputs.
    energy_matrix = coo_matrix(
        (duration[interval] * efficiency[vehicle], (vehicle, entries)),
        shape=(n, peak + 1),
    ).tocsr()
    peak_matrix = coo_matrix(
        (np.r_[np.ones(peak), -np.ones(slots)],
         (np.r_[interval, np.arange(slots)], np.r_[entries, np.full(slots, peak)])),
        shape=(slots, peak + 1),
    ).tocsr()
    limits = [(0, power[i]) for i in vehicle] + [(0, None)]
    ans = linprog(np.r_[np.zeros(peak), 1], A_ub=peak_matrix, b_ub=-load,
                  A_eq=energy_matrix, b_eq=energy, bounds=limits, method="highs",
                  options={"dual_feasibility_tolerance": 1e-9,
                           "primal_feasibility_tolerance": 1e-9})
    if not ans.success:
        raise RuntimeError(f"Minimum charging-only peak, {model['name']}: {ans.message}")
    return float(ans.fun)


def write_csv(path, records):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(records[0]))
        writer.writeheader()
        writer.writerows(records)


def run(replication, output):
    sys.path.insert(0, str(replication))
    from src.energy import optimize

    inputs = sorted((replication / "inputs" / "study").glob("*.json"))
    if len(inputs) != 36:
        raise ValueError(f"Expected 36 study fleets, found {len(inputs)}")
    output.mkdir(parents=True, exist_ok=True)
    multiplier_values = (1.0, 1.1, 1.25, 1.5, None)
    scenarios = {value: [] for value in multiplier_values}
    for path in inputs:
        model = json.loads(path.read_text(encoding="utf-8"))
        if not np.allclose(model["Ulo"], 1) or not np.allclose(model["Uhi"], 1):
            raise ValueError("These experiments require the unit net-import coefficients")
        minimum = minimum_unidirectional_peak(model)
        last_costs = (float("inf"), float("inf"))
        for multiplier in multiplier_values:
            limit = None if multiplier is None else multiplier * minimum
            charge = optimize(model, "charge_only", site_cap=limit)
            robust = optimize(model, "robust_first", site_cap=limit)
            charge_cost, robust_cost = charge["nominal_usd"], robust["worst_usd"]
            saving = charge_cost - robust_cost
            if min(last_costs[0]-charge_cost, last_costs[1]-robust_cost, saving) < -1e-6:
                raise ArithmeticError(f"Cost ordering: {model['name']}, {multiplier}")
            peaks = (charge["protected_peak_kw"], robust["protected_peak_kw"])
            if limit is not None and max(peaks) > limit + 1e-6:
                raise ArithmeticError(f"Physical cap: {model['name']}, {multiplier}")
            last_costs = charge_cost, robust_cost
            scenarios[multiplier].append({
                "case": model["name"], "date": model["metadata"]["horizon_date"],
                "N": model["N"], "q": 1,
                "multiplier": "unbounded" if multiplier is None else f"{multiplier:.2f}",
                "minimum_charging_only_peak_kw": minimum, "limit_kw": limit,
                "charging_only_usd": charge_cost, "v2g_worst_usd": robust_cost,
                "robust_saving_usd": saving, "charging_only_peak_kw": peaks[0],
                "v2g_peak_kw": peaks[1],
                "power_residual_kw": max(charge["power_residual_kw"], robust["power_residual_kw"]),
                "energy_residual_kwh": max(charge["energy_residual_kwh"], robust["energy_residual_kwh"]),
            })
    summary = []
    for rows in scenarios.values():
        charge_total = sum(r["charging_only_usd"] for r in rows)
        robust_total = sum(r["v2g_worst_usd"] for r in rows)
        saving = charge_total - robust_total
        summary.append({
            "multiplier": rows[0]["multiplier"], "fleets": len(rows),
            "charge_only_total_usd": charge_total, "v2g_worst_total_usd": robust_total,
            "total_saving_usd": saving, "pct_total_charging_only": 100*saving/charge_total,
            "positive_fleets": sum(r["robust_saving_usd"] > 1e-6 for r in rows),
        })
    write_csv(output / "common_cap_by_fleet.csv", [r for rows in scenarios.values() for r in rows])
    write_csv(output / "common_cap_summary.csv", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replication", type=Path, default=Path(__file__).resolve().parent,
                        help="Replication root containing inputs/study and src/energy.py")
    parser.add_argument("--out", type=Path, help="Output directory; defaults to replication/results")
    args = parser.parse_args()
    root = args.replication.resolve()
    run(root, args.out or (root / "results"))
