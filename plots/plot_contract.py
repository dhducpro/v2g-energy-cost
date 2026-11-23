"""Draw the one-car example of a discharge withdrawal used in Figure 1."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle


INK = "#192E3A"
MUTED = "#465A63"
PALE = "#C4CFD2"
HIGH = "#9A4123"
LOW = "#1E6272"


def energy_source(ax, x, y, kind):
    """Draw a labeled endpoint; coordinates are expressed in printer's points."""
    def stroke(points, width=1.0):
        ax.plot(*zip(*[(x + a, y + b) for a, b in points]),
                color=INK, lw=width, solid_capstyle="round",
                solid_joinstyle="round", clip_on=False)

    def shape(points):
        ax.add_patch(Polygon([(x + a, y + b) for a, b in points],
                             closed=True, facecolor="white", edgecolor=INK,
                             lw=1, joinstyle="round"))

    if kind == "EV":
        shape([(-10, -1.6), (-9.2, 1.5), (-6.3, 2.1), (-3.7, 5.4),
               (2.5, 5.4), (6.1, 2.1), (8.5, 1.3), (9.7, -1.0),
               (9.7, -3.1), (-10, -3.1)])
        stroke([(-4.8, 2.0), (4.5, 2.0)], .55)
        for position in (-6.0, 5.8):
            ax.add_patch(Circle((x + position, y - 3.1), 1.8,
                                facecolor=INK, edgecolor="white", lw=.35))
        label = "EV"
    elif kind == "Grid":
        stroke([(-5.4, -5.1), (0, 6.0), (5.4, -5.1)])
        stroke([(-6.5, 2.4), (6.5, 2.4)])
        stroke([(-7.4, -.7), (7.4, -.7)])
        stroke([(-3.1, -.8), (3.8, -3.8), (-4.6, -3.8)], .65)
        label = "Grid"
    else:
        shape([(-6.7, -5), (-6.7, 4.5), (4.8, 4.5),
               (4.8, -5)])
        for col in (-4.8, -1.7, 1.4):
            for row in (.7, -2.5):
                ax.add_patch(Rectangle((x + col, y + row), 1.1, 1.4,
                                       color=MUTED, lw=0))
        stroke([(-8.0, -5), (7, -5)], .7)
        label = "Site load"
    ax.text(x, y - 10.8, label, ha="center", va="center",
            fontsize=7.6, color=INK)


def render(output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "ps.fonttype": 42})
    width, height = 238, 166
    fig = plt.figure(figsize=(width / 72, height / 72))
    ax = fig.add_axes((0, 0, 1, 1), xlim=(0, width), ylim=(0, height),
                      aspect="equal")
    ax.set_axis_off()

    def period(y, name, price, color):
        ax.text(5, y, name, fontsize=9.3, fontweight="bold", color=INK,
                va="center")
        ax.text(230, y, price, fontsize=8.3, fontweight="bold", color=color,
                va="center", ha="right")
        ax.plot([5, 230], [y - 9, y - 9], lw=.65, color=PALE)

    def transfer(y, case, source, destination, kwh, color):
        ax.text(5, y - 2, case, va="center", color=INK, fontsize=8.1)
        energy_source(ax, 111, y, source)
        energy_source(ax, 212, y, destination)
        ax.add_patch(FancyArrowPatch((126, y), (197, y), arrowstyle="-|>",
                                     shrinkA=0, shrinkB=0, mutation_scale=8.5,
                                     lw=1.1, color=color))
        ax.text(161, y + 7.2, f"{kwh} kWh", ha="center", va="center",
                color=color, fontsize=7.8)

    period(155, "High-price hour", "30¢/kWh", HIGH)
    transfer(129, "Discharge as planned", "EV", "Site", "1", HIGH)
    transfer(93, "Discharge canceled", "Grid", "Site", "1", HIGH)
    period(58, "Low-price hour", "10¢/kWh", LOW)
    transfer(28, "In both cases", "Grid", "EV", "1.25", LOW)
    fig.savefig(output, metadata={"Title": "How a discharge withdrawal changes a station purchase"})
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("contract_meter.pdf"))
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    render(args.output)
