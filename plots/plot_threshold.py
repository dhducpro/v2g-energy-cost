"""Plot the observed first-zero withdrawal budgets in the fleet study."""
import argparse
import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def load_counts(source):
    counts = Counter()
    with source.open(encoding="utf-8", newline="") as stream:
        cases = set()
        for row in csv.DictReader(stream):
            case, size = row["case"], int(row["N"])
            threshold = int(row["first_zero_observed"])
            if case in cases or not 0 <= threshold <= size:
                raise ValueError(f"Invalid fleet record for {case}")
            cases.add(case)
            counts[size, threshold] += 1
    if not counts:
        raise ValueError("No fleet records were found")
    return counts


def draw(counts, output):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 7.7,
        "pdf.fonttype": 42,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.65,
        "hatch.linewidth": 0.50,
    })
    styles = {
        6: ("#314F6E", "white", ""),
        12: ("#68A8A3", "#335A62", "////"),
        24: ("#DDE4EA", "#778B9A", ""),
    }
    budgets = list(range(max(q for _, q in counts) + 1))
    totals = [sum(counts[size, q] for size in styles) for q in budgets]
    if any(n not in styles for n, _ in counts):
        raise ValueError("Unrecognized fleet size")

    fig, ax = plt.subplots(figsize=(3.45, 2.05))
    fig.subplots_adjust(left=0.145, right=0.985, top=0.835, bottom=0.24)
    bottom = [0] * len(budgets)
    for size, (fill, edge, hatch) in styles.items():
        values = [counts[size, q] for q in budgets]
        present = [k for k, value in enumerate(values) if value]
        ax.bar([budgets[k] for k in present], [values[k] for k in present],
               bottom=[bottom[k] for k in present], width=0.65,
               color=fill, edgecolor=edge, hatch=hatch, linewidth=0.52,
               label=f"{size} vehicles", zorder=2)
        bottom = [a + b for a, b in zip(bottom, values)]
    for q, value in zip(budgets, totals):
        ax.text(q, value + 0.16, str(value), ha="center", va="bottom",
                size=7.2, color="#1D2D39")

    ax.set(xlabel="First withdrawal budget with zero saving",
           ylabel="Number of fleets", xticks=budgets,
           xlim=(budgets[0] - 0.55, budgets[-1] + 0.55),
           ylim=(0, max(totals) + 1.15))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
    ax.set_axisbelow(True)
    ax.grid(axis="y", linewidth=0.45, color="#DCE2E6")
    ax.tick_params(axis="both", length=2.5, width=0.6, pad=3)
    ax.xaxis.labelpad = 6
    ax.yaxis.labelpad = 5
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#70808B")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.015), ncol=3,
              frameon=False, handlelength=1.15, handleheight=1.0,
              borderaxespad=0, columnspacing=1.25, fontsize=7.1)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"Title": "First withdrawal budget with zero guaranteed energy-cost saving"})
    plt.close(fig)


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=here.parent / "results/fleet_thresholds.csv")
    parser.add_argument("--output", type=Path, default=here.parent / "figures/threshold_profile.pdf")
    args = parser.parse_args()
    draw(load_counts(args.input), args.output)
