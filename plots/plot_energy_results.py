"""Reproduce the four empirical figures in the SEGAN energy-cost manuscript."""

import argparse
import csv
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_june import draw_case
from plot_threshold import draw, load_counts

NAVY = "#184C67"
RUST = "#AB4B27"
GRID = "#D6DFE2"
INK = "#22353C"
GREY = "#879DA4"
TOL = 1e-6


def read_rows(filename):
    with filename.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def figure_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 7.4, "pdf.fonttype": 42,
        "axes.labelsize": 7.7, "axes.labelcolor": INK, "text.color": INK,
        "axes.linewidth": 0.55, "axes.edgecolor": "#69808A",
        "axes.spines.right": False, "axes.spines.top": False,
        "xtick.color": INK, "ytick.color": INK,
        "xtick.labelsize": 7.1, "ytick.labelsize": 7.1,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.major.width": 0.55, "ytick.major.width": 0.55,
        "axes.axisbelow": True, "legend.fontsize": 7.1,
    })


def grid(ax):
    ax.grid(axis="y", color=GRID, linewidth=0.46)


def save(fig, output, title):
    fig.savefig(output, metadata={"Title": title})
    plt.close(fig)


def plot_paired(results, output):
    rows = read_rows(results / "daily_monetary_by_fleet.csv")
    if len(rows) != 36 or len({r["case"] for r in rows}) != 36:
        raise ValueError("Expected 36 distinct study fleet-days.")
    pairs = []
    for row in rows:
        base = float(row["charge_only_nominal_usd"])
        nominal = base - float(row["nominal_first_worst_usd"])
        robust = base - float(row["robust_first_worst_usd"])
        pairs.append((0.0 if abs(nominal) <= TOL else nominal,
                      0.0 if abs(robust) <= TOL else robust))
    pairs = sorted(pairs)
    nominal, robust = np.array(pairs).T
    if np.any(robust < nominal - TOL) or np.min(robust) < -TOL:
        raise ValueError("Inconsistent dominance or physical optimum.")
    ranks = np.arange(1, len(pairs) + 1)
    fig, ax = plt.subplots(figsize=(3.45, 2.18))
    fig.subplots_adjust(left=.144, right=.975, bottom=.228, top=.845)
    grid(ax)
    ax.axhline(0, lw=.73, color=INK, zorder=1)
    ax.vlines(ranks, nominal, robust, lw=.8, color=GREY, zorder=2)
    ax.scatter(ranks, nominal, s=11.3, marker="x", lw=.85, color=RUST,
               label="Nominal schedule", zorder=4)
    ax.scatter(ranks, robust, s=11.5, marker="o", lw=.78, color=NAVY,
               facecolor="none", label="Protected schedule", zorder=3)
    ax.set(xlim=(0,37), ylim=(-1.65,5.8), xticks=[1,7,14,21,28,36],
           yticks=[-1,0,1,2,3,4,5], ylabel="Worst-case saving (USD)",
           xlabel="Fleet-day, ordered by nominal-schedule saving")
    ax.legend(loc="lower center", bbox_to_anchor=(.5,1.015), ncol=2,
              frameon=False, handlelength=.8, columnspacing=1.3, borderaxespad=0)
    save(fig, output, "Paired savings with one withdrawal across 36 fleet-days")


def plot_july(results, output):
    rows = sorted(read_rows(results / "july_by_day.csv"), key=lambda r: r["date"])
    dates = [date.fromisoformat(row["date"]) for row in rows]
    days = np.array([d.day for d in dates])
    if list(days) != list(range(1, 32)) or any(d.year != 2019 or d.month != 7 for d in dates):
        raise ValueError("Expected every calendar day of July 2019 in order.")
    values = np.array([float(r["charge_only_nominal_usd"])-float(r["robust_first_worst_usd"])
                       for r in rows])
    values[np.abs(values) <= TOL] = 0
    if values.min() < -TOL:
        raise ValueError("Protected savings must be nonnegative.")
    weekend = np.array([d.weekday() >= 5 for d in dates])
    fig, ax = plt.subplots(figsize=(3.45,2.04))
    fig.subplots_adjust(left=.147, right=.978, bottom=.235, top=.846)
    grid(ax)
    # Mark zero outcomes explicitly, including the weekday that ties the baseline.
    for sel, color, marker, label in [(~weekend, NAVY, "o", "Weekday"),
                                      (weekend, RUST, "s", "Weekend")]:
        ax.vlines(days[sel], 0, values[sel], color=color, lw=1.0)
        ax.scatter(days[sel], values[sel], s=13.0, color=color, lw=.75,
                   facecolor="white" if marker=="s" else color, marker=marker,
                   zorder=4, label=label)
    ax.axhline(0,lw=.65,color=INK)
    ax.set(xlim=(.1,31.9), ylim=(-.75,10.1), xticks=[1,5,10,15,20,25,31],
           yticks=[0,2,4,6,8,10], xlabel="Day of July 2019",
           ylabel="Guaranteed saving (USD)")
    ax.legend(loc="lower center", bbox_to_anchor=(.5,1.012), frameon=False,
              ncol=2, handlelength=.8, columnspacing=1.5, borderaxespad=0)
    save(fig, output, "Daily guaranteed energy saving across July 2019")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="Path to the replication folder with inputs/, results/ and src/.")
    parser.add_argument("--output-dir", type=Path, help="Directory for the four empirical PDFs.")
    parser.add_argument("--results", type=Path, help="Results folder; defaults to results/ at the replication root.")
    args = parser.parse_args()
    results = args.results or args.root / "results"
    output = args.output_dir or args.root / "figures"
    output.mkdir(parents=True, exist_ok=True)
    draw_case(args.root, output / "june_example.pdf")
    draw(load_counts(results / "fleet_thresholds.csv"), output / "threshold_profile.pdf")
    figure_style()
    plot_paired(results, output / "fleet_savings_pairs.pdf")
    plot_july(results, output / "july_daily_savings.pdf")


if __name__ == "__main__":
    main()
