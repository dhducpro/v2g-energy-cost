from pathlib import Path
import csv, argparse
def compute(records, unit_charge):
    fixed, selected, count = 0., 0., 0
    for r in records:
        saving = float(r["charge_only_nominal_usd"])-float(r["robust_first_worst_usd"])
        allowance = unit_charge*float(r["robust_first_discharge_kwh"])
        margin = saving-allowance
        fixed += margin
        selected += max(0,margin)
        count += margin > 1e-6
    return fixed, selected, count
if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Ex-post sensitivity of the energy-optimal schedules to an illustrative discharge payment.")
    parser.add_argument("results",type=Path,help="Directory with results/daily_monetary_by_fleet.csv or the CSV directory itself.")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    f=args.results/"daily_monetary_by_fleet.csv"
    with f.open(encoding="utf8",newline="") as stream:
        rows=list(csv.DictReader(stream))
    assert len(rows)==36
    with args.out.open("w",encoding="utf8",newline="") as stream:
        dst=csv.writer(stream)
        dst.writerow(["usd_per_scheduled_discharge_kwh","fixed_schedules_margin_usd","choose_each_day_lower_bound_usd","strictly_positive_fleets"])
        for unit in [0.,0.01,0.03,0.05,0.1]:
            fixed,chosen,count=compute(rows,unit)
            dst.writerow([unit,fixed,chosen,count])
            print(f"{unit:.2f}   {fixed:.3f}   {chosen:.3f}   {count}")
