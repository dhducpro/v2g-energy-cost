"""Optimize the paper's energy-only schedules on its supplied numerical inputs."""
import argparse
import csv
from datetime import date, timedelta
import json
from pathlib import Path

import numpy as np
from src.energy import optimize, prices_for_day


ROOT = Path(__file__).resolve().parent


def write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def cases(directory):
    for path in sorted(directory.glob("*.json")):
        yield json.loads(path.read_text(encoding="utf-8"))


def daily_record(m, full):
    weighted_prices = prices_for_day(m) * np.asarray(m["duration"])
    background_cost = np.dot(weighted_prices, m["L"])
    record = {
        "case": m["name"],
        "date": m["metadata"]["horizon_date"],
        "N": m["N"],
        "q": 1,
        "background_energy_cost_usd": background_cost,
    }
    modes = ["charge_only", "robust_first"]
    if full:
        modes = ["charge_only", "nominal_first", "robust_first"]
    for mode in modes:
        for metric, value in optimize(m, mode).items():
            record[f"{mode}_{metric}"] = value

    charging_only_cost = record["charge_only_nominal_usd"]
    for mode in modes:
        for bill in ("nominal_usd", "worst_usd"):
            key = f"{mode}_{bill}"
            record[key + "_saving"] = charging_only_cost - record[key]
            record[key + "_net_of_background"] = record[key] - background_cost
    return record


def study(output):
    rows = [daily_record(model, full=True) for model in cases(ROOT / "inputs/study")]
    assert len(rows) == 36 and sum(row["N"] for row in rows) == 504
    write(output / "daily_monetary_by_fleet.csv", rows)


def seasons(output):
    rows = []
    tariff_weeks = [("winter", date(2019, 4, 1)), ("summer", date(2019, 7, 1))]
    for season, monday in tariff_weeks:
        for model in cases(ROOT / "inputs/study"):
            actual_date = model["metadata"]["horizon_date"]
            # Keep the fleet's day of the week when switching tariff seasons.
            tariff_date = monday + timedelta(days=date.fromisoformat(actual_date).weekday())
            prices = prices_for_day({**model, "metadata": {"horizon_date": tariff_date.isoformat()}})
            charging_only = optimize(model, "charge_only", price_vector=prices)["worst_usd"]
            v2g = optimize(model, "robust_first", price_vector=prices)
            saving = charging_only - v2g["worst_usd"]
            background_energy = np.asarray(model["L"]) * model["duration"]
            rows.append({
                "season": season,
                "case": model["name"],
                "date": actual_date,
                "N": model["N"],
                "charge_only_usd": charging_only,
                "v2g_worst_usd": v2g["worst_usd"],
                "saving_usd": saving,
                "saving_pct": 100 * saving / charging_only,
                "background_usd": np.dot(background_energy, prices),
                "v2g_discharge_kwh": v2g["discharge_kwh"],
            })
    assert len(rows) == 72
    write(output / "same_fleets_by_season.csv", rows)


def july(output):
    rows = [daily_record(model, full=False) for model in cases(ROOT / "inputs/july")]
    sessions = sum(row["N"] for row in rows)
    assert len(rows) == 31 and sessions == 578
    write(output / "july_by_day.csv", rows)
    charging_only = sum(row["charge_only_nominal_usd"] for row in rows)
    v2g = sum(row["robust_first_worst_usd"] for row in rows)
    saving = charging_only - v2g
    totals = {
        "days": len(rows),
        "sessions": sessions,
        "charging_bill_usd": charging_only,
        "v2g_bill_usd": v2g,
        "savings_usd": saving,
        "savings_percent": 100 * saving / charging_only,
    }
    write(output / "july_totals.csv", [totals])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", choices=("study", "seasons", "july", "all"),
                        default="all", nargs="?")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    for name, run in [("study", study), ("seasons", seasons), ("july", july)]:
        if args.experiment in (name, "all"):
            run(args.out)
            print(f"Completed {name}.", flush=True)
