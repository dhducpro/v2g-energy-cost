"""Rebuild the June 5 vehicle-withdrawal meter profiles from the archived numerical input."""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


def reconstruct(model, mode, solver, audit):
    schedule = {}

    def record_schedule(instance, q, charge, discharge, peak):
        schedule["charge"] = charge.copy()
        schedule["discharge"] = discharge.copy()
        return audit(instance, q, charge, discharge, peak)

    result = solver.optimize(model, mode, audit=record_schedule)
    duration = np.array(model["duration"])
    prices = solver.prices_for_day(model)
    discharge = schedule["discharge"]
    net = np.array(model["L"]) + schedule["charge"].sum(axis=0) - discharge.sum(axis=0)
    delivered = np.maximum(0, net)
    withdrawals = np.maximum(0, net[None, :] + discharge)
    bills = withdrawals @ (duration * prices)
    vehicle = 1  # The same second vehicle in each policy, indexed from zero.
    nominal_bill = delivered @ (duration * prices)
    missing_kwh = (withdrawals[vehicle] - delivered) @ duration
    assert abs(bills[vehicle] - bills.max()) < 1e-6
    assert abs(bills[vehicle] - result["worst_usd"]) < 1e-6
    assert abs(nominal_bill - result["nominal_usd"]) < 1e-6
    assert abs(bills[vehicle] - nominal_bill - 0.267 * missing_kwh) < 1e-6
    # Confirm the printed case figures match the reconstructed meter profile.
    expected = (31.59, 35.12, 13.20) if mode == "nominal_first" else (33.04, 34.05, 3.80)
    assert tuple(round(v, 2) for v in (nominal_bill, bills[vehicle], missing_kwh)) == expected
    return delivered, withdrawals[vehicle], missing_kwh


def plot_panel(ax, edges, delivered, withdrawn, title, missing_kwh):
    navy, rust, pale = "#164B68", "#A33D20", "#F3EDDA"
    ax.axvspan(12, 18, color=pale, lw=0, zorder=0)
    ax.stairs(delivered, edges, baseline=None, color=navy, lw=1.55, zorder=3)
    # Show dashed segments where the scenarios separate. A single solid trace
    # represents both cases wherever meter import coincides.
    changed = withdrawn - delivered > 1e-6
    starts = np.flatnonzero(np.diff(np.r_[False, changed, False].astype(int)) == 1)
    ends = np.flatnonzero(np.diff(np.r_[False, changed, False].astype(int)) == -1)
    for start, end in zip(starts, ends):
        heights = withdrawn[start:end]
        x = np.repeat(edges[start:end + 1], 2)[1:-1]
        y = np.repeat(heights, 2)
        ax.fill_between(x, np.repeat(delivered[start:end], 2), y,
                        color="#E7A88A", alpha=0.49, lw=0, zorder=1)
        before = withdrawn[start - 1] if start else delivered[0]
        after = withdrawn[end] if end < len(withdrawn) else delivered[-1]
        ax.plot(np.r_[x[0], x, x[-1]], np.r_[before, y, after],
                color=rust, lw=1.7, dashes=(3.2, 2.0), zorder=4)
    ax.set_xlim(8, 19)
    ax.set_ylim(-0.6, 30)
    ax.set_yticks([0, 10, 20, 30])
    ax.set_ylabel("Import (kW)", fontsize=7.6, labelpad=2.0)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#B5BABD", linewidth=0.45, alpha=0.65)
    ax.spines[["top", "right"]].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color("#79838A")
        ax.spines[spine].set_linewidth(0.65)
    ax.tick_params(axis="both", labelsize=7.1, color="#79838A", width=0.6, length=3)
    ax.set_title(title, loc="left", fontsize=7.9, fontweight="bold", pad=4)
    ax.set_title(f"{missing_kwh:.2f} kWh replaced", loc="right", fontsize=7.3, pad=4)


def draw_case(root, output):
    sys.path.insert(0, str(root/"src"))
    import energy as solver
    from energy import audit_schedule
    path = root/"inputs"/"study"/"A201906_04_N06.json"
    model = json.loads(path.read_text(encoding="utf-8"))
    assert model["metadata"]["horizon_date"] == "2019-06-05" and model["N"] == 6
    edges = np.asarray(model["boundaries_seconds"]) / 3600
    nominal = reconstruct(model, "nominal_first", solver, audit_schedule)
    robust = reconstruct(model, "robust_first", solver, audit_schedule)
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "axes.labelcolor": "#26343C", "text.color": "#26343C",
                         "xtick.color": "#26343C", "ytick.color": "#26343C"})
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(3.45, 2.75))
    for ax, title, data in zip(axes, ["(a)  Nominal V2G", "(b)  Robust V2G"], [nominal, robust]):
        plot_panel(ax, edges, data[0], data[1], title, data[2])
    axes[1].set_xticks([8, 10, 12, 14, 16, 18],
                       ["08:00","10:00","12:00","14:00","16:00","18:00"])
    axes[1].set_xlabel("Local time", fontsize=7.6, labelpad=2)
    legend = [Line2D([], [], color="#164B68", lw=1.55, label="All deliver"),
              Line2D([], [], color="#A33D20", lw=1.7, dashes=(3.2,2.0), label="Vehicle 2 withdraws"),
              Patch(facecolor="#F3EDDA", edgecolor="#C9BC94", linewidth=0.45, label="On-peak")]
    fig.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.52,1.018), ncol=3,
               frameon=False, fontsize=7.1, handlelength=1.6, columnspacing=0.9)
    fig.subplots_adjust(left=0.154, right=0.976, bottom=0.137, top=0.84, hspace=0.54)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"Title":"June 5 six-EV station meter schedules"})
    plt.close(fig)


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=root)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    draw_case(args.root, args.output or args.root / "figures" / "june_example.pdf")
