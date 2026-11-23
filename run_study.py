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
    with path.open("w",newline="",encoding="utf-8") as output:
        writer = csv.DictWriter(output,fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def cases(directory):
    for path in sorted(directory.glob("*.json")):
        yield json.loads(path.read_text(encoding="utf-8"))


def daily_record(m, full):
    record = dict(case=m["name"], date=m["metadata"]["horizon_date"], N=m["N"], q=1)
    record["background_energy_cost_usd"] = np.dot(prices_for_day(m)*np.array(m["duration"]),m["L"])
    modes = ["charge_only","nominal_first","robust_first"] if full else ["charge_only","robust_first"]
    for mode in modes:
        outcomes = optimize(m, mode)
        record.update({f"{mode}_{k}":v for k,v in outcomes.items()})
    for mode in modes:
        for label in ["nominal_usd","worst_usd"]:
            key = f"{mode}_{label}"
            record[key+"_saving"] = record["charge_only_nominal_usd"]-record[key]
            record[key+"_net_of_background"] = record[key]-record["background_energy_cost_usd"]
    return record


def study(output):
    rows = [daily_record(m,True) for m in cases(ROOT/"inputs/study")]
    assert len(rows)==36 and sum(r["N"] for r in rows)==504
    write(output/"daily_monetary_by_fleet.csv",rows)


def seasons(output):
    rows = []
    for season,monday in [("winter",date(2019,4,1)),("summer",date(2019,7,1))]:
        for m in cases(ROOT/"inputs/study"):
            actual = m["metadata"]["horizon_date"]
            date_for_rates = monday + timedelta(days=date.fromisoformat(actual).weekday())
            prices = prices_for_day({**m,"metadata":{"horizon_date":date_for_rates.isoformat()}})
            baseline = optimize(m,"charge_only",price_vector=prices)
            v2g = optimize(m,"robust_first",price_vector=prices)
            saving = baseline["worst_usd"] - v2g["worst_usd"]
            rows.append(dict(season=season,case=m["name"],date=actual,N=m["N"],
                             charge_only_usd=baseline["worst_usd"],v2g_worst_usd=v2g["worst_usd"],
                             saving_usd=saving,saving_pct=100*saving/baseline["worst_usd"],
                             background_usd=np.dot(np.array(m["L"])*m["duration"],prices),
                             v2g_discharge_kwh=v2g["discharge_kwh"]))
    assert len(rows)==72
    write(output/"same_fleets_by_season.csv",rows)


def july(output):
    rows = [daily_record(m,False) for m in cases(ROOT/"inputs/july")]
    assert len(rows)==31 and sum(r["N"] for r in rows)==578
    write(output/"july_by_day.csv",rows)
    charge = sum(r["charge_only_nominal_usd"] for r in rows)
    worst = sum(r["robust_first_worst_usd"] for r in rows)
    write(output/"july_totals.csv",[dict(days=len(rows),sessions=sum(r["N"] for r in rows),
          charging_bill_usd=charge,v2g_bill_usd=worst,savings_usd=charge-worst,
          savings_percent=100*(charge-worst)/charge)])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment",choices=("study","seasons","july","all"),default="all",nargs="?")
    parser.add_argument("--out",type=Path,default=ROOT/"results")
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    for name,func in [("study",study),("seasons",seasons),("july",july)]:
        if args.experiment in (name,"all"):
            func(args.out)
            print(f"Completed {name}.",flush=True)
