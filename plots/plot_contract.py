"""Render the schematic for the withdrawal contract and the two energy bills."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch


def render(output):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 7.2,
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    ink = "#182B39"
    muted = "#4A5961"
    charge = "#14678A"
    discharge = "#AB5427"
    neutral = "#CBD2D4"
    fig = plt.figure(figsize=(3.36, 2.42))
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()

    def write(x, y, value, *, size=7.2, color=ink, ha="left", bold=False):
        return ax.text(x, y, value, ha=ha, va="center", fontsize=size,
                       color=color, fontweight="bold" if bold else "normal",
                       linespacing=1.20)

    def rule(x0, x1, y, color=neutral, width=0.65):
        ax.plot((x0, x1), (y, y), color=color, lw=width, solid_capstyle="butt")

    def flow(x, y, start, end, flow_color):
        write(x-0.13, y, start, color=flow_color)
        write(x+0.13, y, end, color=flow_color, ha="right")
        arrow = FancyArrowPatch((x-0.027, y), (x+0.039, y), arrowstyle="-|>",
                                color=flow_color, lw=1.0, mutation_scale=7.5)
        ax.add_patch(arrow)

    write(0.015, 0.954, "(a)  Power delivered during one session", size=8.0, bold=True)
    write(0.485, 0.817, "Discharge interval", size=7.1, ha="center")
    write(0.819, 0.817, "Recharge interval", size=7.1, ha="center")
    rule(0.018, 0.978, 0.771)
    write(0.024, 0.675, "Permission kept")
    write(0.024, 0.561, "Withdrawn")
    flow(0.485, 0.675, "EV", "Site", discharge)
    write(0.485, 0.561, "No discharge", color=muted, ha="center")
    flow(0.819, 0.675, "Site", "EV", charge)
    flow(0.819, 0.561, "Site", "EV", charge)
    rule(0.018, 0.978, 0.49)

    write(0.015, 0.394, "(b)  Station bill under the same tariff", size=8.0, bold=True)
    write(0.020, 0.277, "Optimized charging only")
    write(0.020, 0.153, "Protected V2G schedule")
    write(0.563, 0.277, r"$C_U$", size=8.0, ha="center", bold=True)
    write(0.563, 0.153, r"$W_q$", size=8.0, ha="center", bold=True)
    write(0.637, 0.277, "Minimum bill", size=7.0)
    write(0.637, 0.153, "Min. worst-case bill", size=7.0)
    write(0.020, 0.043, "Each bill prices positive net imports interval by interval.",
          size=7.0, color=muted)
    fig.savefig(output, metadata={"Title": "Discharge withdrawal and station energy bills"})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("contract_meter.pdf"))
    args = parser.parse_args()
    args.output.parent.mkdir(exist_ok=True, parents=True)
    render(args.output)
