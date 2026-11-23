"""Two-hour surplus test for the revised model with zero export credit."""
import argparse
import csv
from pathlib import Path

from mode_check import control


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("surplus_mode_check.csv"))
    args=parser.parse_args()
    model=dict(name="two_hour_surplus", N=2, T=2, availability=[[1,1],[1,1]],
               duration=[1.,1.], L=[-4.,3.], Pc=[1.,1.], Pd=[1.,1.], eta_c=[.9,.9],
               eta_d=[.9,.9], E0=[2.,2.], Etarget=[2.,2.], Emin=[0.,0.], Emax=[2.9,2.9])
    prices=[.10,.30]
    only=control({**model, "Pd":[0.,0.]},False,20,prices)
    lp=control(model,False,20,prices)
    mip=control(model,True,20,prices)
    for sol, expected in [(only,.9),(lp,.657),(mip,.657)]:
        assert sol["status"]==0 and abs(sol["worst_bill_usd"]-expected)<1e-8
    assert abs(mip["delivered_min_net_kw"]+2)<1e-8
    row=dict(charge_only_usd=only["worst_bill_usd"],lp_usd=lp["worst_bill_usd"],
             binary_usd=mip["worst_bill_usd"],binary_minus_lp_usd=mip["objective_usd"]-lp["objective_usd"],
             saving_usd=only["worst_bill_usd"]-mip["worst_bill_usd"],
             binary_full_delivery_first_net_kw=mip["delivered_min_net_kw"],
             binary_full_delivery_export_kwh=mip["delivered_export_kwh"],
             binary_battery_violation_kwh=mip["battery_violation_kwh"],binary_overlap_kw=mip["overlap_kw"])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=row.keys())
        writer.writeheader()
        writer.writerow(row)
    print(row)


if __name__=="__main__":
    main()
