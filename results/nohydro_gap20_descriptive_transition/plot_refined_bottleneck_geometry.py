#!/usr/bin/env python3
"""Plot the exact production boundary, exceptional links, and indexed crossings."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyze_approach_current_balance import load_gap20_production  # noqa: E402
from run_bottleneck_commitment_branches import SOURCE  # noqa: E402


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def plot_boundary(ax, mod, size: float = 0.6) -> None:
    ax.scatter(mod.x_coord, mod.y_coord, s=size, color="black", alpha=0.75,
               rasterized=True, label="production boundary")
    ax.set_aspect("equal")
    ax.grid(alpha=0.18)


def intersects_arc(a: np.ndarray, b: np.ndarray, center: np.ndarray, radius: float,
                   theta_min: float, theta_max: float) -> bool:
    """Return whether segment a-b intersects the specified production circle arc."""
    direction = b - a; offset = a - center
    roots = np.roots((direction @ direction, 2 * offset @ direction,
                      offset @ offset - radius * radius))
    for root in roots:
        if abs(root.imag) > 1e-12 or not -1e-12 <= root.real <= 1 + 1e-12:
            continue
        point = a + root.real * direction
        theta = np.arctan2(point[1] - center[1], point[0] - center[0])
        if theta_min - 1e-10 <= theta <= theta_max + 1e-10:
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--input", type=Path,
                        default=HERE / "outputs" / "bottleneck_directional_refined")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output or args.input / "geometry_visibility_review.png"
    mod = load_gap20_production(args.source.resolve())
    links = read_csv(args.input / "offcentral_link_endpoints.csv")
    points = read_csv(args.input / "offcentral_snapshots.csv")
    indexed = read_csv(args.input / "indexed_first_crossing_audit.csv")

    unique_links: dict[tuple[str, int, int], dict[str, str]] = {}
    for row in links:
        key = (row["sample_id"], *sorted((int(row["fish_i"]), int(row["fish_j"]))))
        unique_links.setdefault(key, row)
    for row in unique_links.values():
        a = np.array([float(row["xi"]), float(row["yi"])])
        b = np.array([float(row["xj"]), float(row["yj"])])
        left = intersects_arc(a, b, mod.circle_center2, mod.delta, 0, mod.th_1)
        right = intersects_arc(a, b, mod.circle_center1, mod.delta, mod.th_8, np.pi)
        if not (left and right):
            raise RuntimeError("Exceptional sightline does not intersect both upper boundary arcs")

    fig = plt.figure(figsize=(13.2, 9.2))
    grid = fig.add_gridspec(2, 2, height_ratios=(1.08, 1), hspace=0.28, wspace=0.22)
    ax_full = fig.add_subplot(grid[0, 0])
    ax_hist = fig.add_subplot(grid[0, 1])
    ax_left = fig.add_subplot(grid[1, 0])
    ax_right = fig.add_subplot(grid[1, 1])

    plot_boundary(ax_full, mod, 0.45)
    ax_full.axhspan(-mod.gap / 2, mod.gap / 2, color="C2", alpha=0.10,
                    label="nominal central aperture")
    ax_full.set(xlim=(-7.5, 7.5), ylim=(-7.5, 7.5), xlabel="x", ylabel="y",
                title="A  Exact production boundary")
    ax_full.legend(loc="lower right", fontsize=8)

    crossing_y = np.asarray([float(row["max_abs_y"]) for row in indexed
                             if row["max_abs_y"] not in ("", "None")])
    ax_hist.hist(crossing_y, bins=np.linspace(0, mod.gap / 2, 24), color="C2", alpha=0.75)
    ax_hist.axvline(mod.gap / 2, color="black", linestyle="--", linewidth=1,
                    label=f"gap/2 = {mod.gap / 2:.3f}")
    ax_hist.set(xlabel="|y| at indexed first-crossing step", ylabel="indexed episodes",
                title=f"B  All {len(indexed)} indexed first crossings are central")
    ax_hist.grid(alpha=0.18); ax_hist.legend(fontsize=8)

    panels = ((ax_left, "pre_s1_k985000", "C3", "C"),
              (ax_right, "pre_s12_k240000", "C0", "D"))
    for ax, sample, color, letter in panels:
        plot_boundary(ax, mod, 0.8)
        sample_points = [row for row in points if row["sample_id"] == sample]
        ax.scatter([float(row["x"]) for row in sample_points],
                   [float(row["y"]) for row in sample_points], s=12,
                   color="0.68", zorder=2, label="fish")
        seen: set[tuple[int, int]] = set()
        for row in links:
            if row["sample_id"] != sample:
                continue
            pair = tuple(sorted((int(row["fish_i"]), int(row["fish_j"]))))
            if pair in seen:
                continue
            seen.add(pair)
            ax.plot([float(row["xi"]), float(row["xj"])],
                    [float(row["yi"]), float(row["yj"])], color=color,
                    linewidth=1.2, alpha=0.9, zorder=3)
        ax.set(xlim=(-2, 2), ylim=(5.75, 7.25), xlabel="x", ylabel="y",
               title=f"{letter}  {sample}: {len(seen)} exceptional links")
        ax.text(0.03, 0.04, "Sightlines cross the rounded wall arcs",
                transform=ax.transAxes, fontsize=8,
                bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.9))

    fig.suptitle("Refined geometry audit: native off-central visibility is through-wall",
                 fontsize=14, y=0.98)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
