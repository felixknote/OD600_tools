"""
01_plate_heatmap_edge_effect.py
================================
#env = mAIcrobe

OD600 plate heatmap and edge-effect analysis for 96-well plates.

Reads endpoint OD600 measurements from multiple plates (CSV, semicolon-delimited,
comma decimal separator). Each file represents one plate (biological × technical
replicate combination), named P<bio>_<tech>.csv (e.g. P1_1.csv, P2_1.csv).

Output:  <DATA_DIR>/plate_analysis/
    ├── P1_1.png, P1_2.png, …        individual plate heatmaps
    ├── average_plate.png             mean across all plates
    ├── correlation_od_vs_distance.png violin + stats for edge-distance effect
    └── zone_map.png                  96-well layout coloured by edge distance
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from scipy import stats
from scipy.stats import gaussian_kde

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

class Config:
    DATA_DIR   = Path(__file__).parent
    OUTPUT_DIR = DATA_DIR / "plate_analysis"
    PLATE_GLOB = "P*.csv"       # matches P1_1.csv, P2_1.csv, P3_2.csv …
    CMAP       = "Reds"
    DPI        = 300

ROW_LABELS = list("ABCDEFGH")
COL_LABELS = [str(c) for c in range(1, 13)]

# Edge-distance color scheme (0 = outer edge → 3 = center)
DIST_COLORS: dict[int, str] = {
    0: "#c0392b",
    1: "#e07070",
    2: "#a8c8e8",
    3: "#2980b9",
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_plates() -> dict[str, np.ndarray]:
    """Load all plate CSVs into a {stem: (8,12) array} dict."""
    plate_files = sorted(Config.DATA_DIR.glob(Config.PLATE_GLOB))
    if not plate_files:
        raise FileNotFoundError(f"No files matching '{Config.PLATE_GLOB}' in {Config.DATA_DIR}")
    plates = {}
    for f in plate_files:
        df = pd.read_csv(f, sep=";", header=None, decimal=",")
        plates[f.stem] = df.values.astype(float)
    return plates


def build_tidy_df(plates: dict[str, np.ndarray]) -> pd.DataFrame:
    """Build long-form DataFrame with edge-distance annotation for each well."""
    rows_idx, cols_idx = np.meshgrid(np.arange(8), np.arange(12), indexing="ij")
    zone_grid = np.minimum(
        np.minimum(rows_idx, 7 - rows_idx),
        np.minimum(cols_idx, 11 - cols_idx),
    )
    records = []
    for name, data in plates.items():
        bio, tech = name.split("_")
        for r in range(8):
            for c in range(12):
                records.append({
                    "plate":     name,
                    "bio_rep":   bio,
                    "tech_rep":  tech,
                    "row":       r,
                    "col":       c,
                    "edge_dist": int(zone_grid[r, c]),
                    "OD600":     data[r, c],
                })
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _format_plate_axes(ax: plt.Axes) -> None:
    ax.set_xticks(range(12))
    ax.set_xticklabels(COL_LABELS, fontsize=8)
    ax.set_yticks(range(8))
    ax.set_yticklabels(ROW_LABELS, fontsize=8)
    ax.tick_params(length=0)


def plot_plate_heatmap(
    data: np.ndarray,
    title: str,
    save_path: Path,
    vmin: float,
    vmax: float,
) -> None:
    """Save a single 96-well plate heatmap."""
    fig, ax = plt.subplots(figsize=(6, 3.5))
    im = ax.imshow(data, cmap=Config.CMAP, vmin=vmin, vmax=vmax, aspect="auto")
    _format_plate_axes(ax)
    ax.set_title(title, fontsize=11, pad=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("OD600", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(save_path, dpi=Config.DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {save_path.name}")


def _sig_label(p: float) -> str:
    return "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))


def plot_edge_correlation(df_all: pd.DataFrame) -> None:
    """Violin plot of OD600 by edge distance with pairwise Mann-Whitney brackets."""
    distances   = sorted(df_all["edge_dist"].unique())
    plate_means = df_all.groupby(["plate", "edge_dist"])["OD600"].mean().reset_index()
    r_corr, p_pearson = stats.pearsonr(df_all["edge_dist"], df_all["OD600"])
    rng = np.random.default_rng(42)

    # Build KDE-extended violin stats
    vpstats = []
    violin_ymax = -np.inf
    for dist in distances:
        vals = df_all.loc[df_all["edge_dist"] == dist, "OD600"].values
        kde  = gaussian_kde(vals)
        bw   = kde.factor * vals.std()
        coords = np.linspace(vals.min() - 2 * bw, vals.max() + 2 * bw, 300)
        vpstats.append({
            "coords": coords, "vals": kde(coords),
            "mean": vals.mean(), "median": float(np.median(vals)),
            "min": vals.min(), "max": vals.max(),
        })
        violin_ymax = max(violin_ymax, coords.max())

    fig, ax = plt.subplots(figsize=(9, 5))
    vp = ax.violin(vpstats, positions=list(range(len(distances))),
                   widths=0.7, showmedians=False, showextrema=False)
    for body, dist in zip(vp["bodies"], distances):
        body.set_facecolor(DIST_COLORS[dist])
        body.set_alpha(0.75)

    means_y = []
    for pos, dist in enumerate(distances):
        vals = df_all.loc[df_all["edge_dist"] == dist, "OD600"].values
        color = DIST_COLORS[dist]
        ax.scatter(pos + rng.uniform(-0.14, 0.14, size=len(vals)), vals,
                   color=color, s=6, alpha=0.35, linewidths=0, zorder=2)
        pm = plate_means.loc[plate_means["edge_dist"] == dist, "OD600"].values
        ax.scatter(pos + rng.uniform(-0.07, 0.07, size=len(pm)), pm,
                   color="black", s=40, zorder=5, marker="D",
                   linewidths=0.5, edgecolors="white")
        means_y.append(float(vals.mean()))

    ax.plot(range(len(distances)), means_y, color="black", linewidth=1.5,
            zorder=6, marker="o", markersize=5)

    # Significance brackets between adjacent distances
    span      = df_all["OD600"].max() - df_all["OD600"].min()
    y_bracket = violin_ymax + span * 0.04
    tick_h    = span * 0.02
    for i in range(len(distances) - 1):
        a, b = distances[i], distances[i + 1]
        va = df_all.loc[df_all["edge_dist"] == a, "OD600"].values
        vb = df_all.loc[df_all["edge_dist"] == b, "OD600"].values
        _, p_ab = stats.mannwhitneyu(va, vb, alternative="two-sided")
        label   = _sig_label(p_ab)
        color   = "black" if label != "ns" else "#999999"
        ax.plot([a, b], [y_bracket, y_bracket], color=color, linewidth=1, clip_on=False)
        ax.plot([a, a], [y_bracket - tick_h, y_bracket], color=color, linewidth=1, clip_on=False)
        ax.plot([b, b], [y_bracket - tick_h, y_bracket], color=color, linewidth=1, clip_on=False)
        ax.text((a + b) / 2, y_bracket + tick_h * 0.3, label,
                ha="center", va="bottom", fontsize=8, color=color, clip_on=False)

    ax.set_ylim(df_all["OD600"].min() - span * 0.04, y_bracket + span * 0.10)
    ax.set_xticks(range(len(distances)))
    ax.set_xticklabels([str(d) for d in distances], fontsize=10)
    ax.set_xlabel("Edge Distance  (0 = outer edge)", fontsize=10)
    ax.set_ylabel("OD600", fontsize=10)
    ax.set_title("OD600 Edge Distance Correlation", fontsize=10, pad=8)

    legend_handles = [
        Patch(facecolor=DIST_COLORS[d], alpha=0.75,
              label=f"{d}" + (" (edge)" if d == 0 else
                              " (center)" if d == max(distances) else ""))
        for d in distances
    ]
    legend_handles += [
        plt.Line2D([0], [0], color="none", marker="D", markerfacecolor="black",
                   markersize=6, label="Per-plate means"),
        plt.Line2D([0], [0], color="black", linewidth=1.5, marker="o", markersize=5,
                   label=f"r = {r_corr:.3f}  (p = {p_pearson:.2e})"),
    ]
    ax.legend(handles=legend_handles, fontsize=8, framealpha=0.8,
              bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0,
              title="Edge Distance", title_fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    out = Config.OUTPUT_DIR / "correlation_od_vs_distance.png"
    fig.savefig(out, dpi=Config.DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def plot_zone_map() -> None:
    """96-well plate layout coloured by edge distance (0 = edge, 3 = center)."""
    rows_idx, cols_idx = np.meshgrid(np.arange(8), np.arange(12), indexing="ij")
    zone_grid = np.minimum(
        np.minimum(rows_idx, 7 - rows_idx),
        np.minimum(cols_idx, 11 - cols_idx),
    ).astype(float)

    cmap_zones = ListedColormap([DIST_COLORS[d] for d in sorted(DIST_COLORS)])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], len(DIST_COLORS))

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.imshow(zone_grid, cmap=cmap_zones, norm=norm, aspect="auto")
    _format_plate_axes(ax)

    for r in range(8):
        for c in range(12):
            ax.text(c, r, str(int(zone_grid[r, c])),
                    ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")

    legend_handles = [
        Patch(facecolor=DIST_COLORS[d], alpha=0.9,
              label=f"{d}" + (" — outer edge" if d == 0 else
                              " — center"    if d == max(DIST_COLORS) else ""))
        for d in sorted(DIST_COLORS)
    ]
    ax.legend(handles=legend_handles, fontsize=8, framealpha=0.9,
              bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0,
              title="Edge Distance", title_fontsize=8)
    ax.set_title("Plate zone map", fontsize=11, pad=8)
    fig.tight_layout()

    out = Config.OUTPUT_DIR / "zone_map.png"
    fig.savefig(out, dpi=Config.DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out.name}")


def print_summary(df_all: pd.DataFrame) -> None:
    r_corr, p_pearson = stats.pearsonr(df_all["edge_dist"], df_all["OD600"])
    print(f"\n{'Dist':<6}  {'Mean':>7}  {'SD':>7}  {'n':>5}")
    print("-" * 32)
    for dist_val, grp in df_all.groupby("edge_dist"):
        print(f"{dist_val:<6}  {grp['OD600'].mean():>7.4f}  {grp['OD600'].std():>7.4f}  {len(grp):>5}")
    print(f"\nPearson r (OD vs edge distance): {r_corr:.4f}  p = {p_pearson:.3g}")
    max_dist = int(df_all["edge_dist"].max())
    _, p_0max = stats.mannwhitneyu(
        df_all.loc[df_all["edge_dist"] == 0, "OD600"].values,
        df_all.loc[df_all["edge_dist"] == max_dist, "OD600"].values,
        alternative="two-sided",
    )
    print(f"Mann-Whitney U (edge dist 0 vs {max_dist}): p = {p_0max:.3g}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("=== OD600 Plate Heatmap & Edge-Effect Analysis ===")
    Config.OUTPUT_DIR.mkdir(exist_ok=True)

    # Load
    plates = load_plates()
    print(f"  Loaded {len(plates)} plate(s): {', '.join(sorted(plates))}")

    all_values = np.concatenate([v.ravel() for v in plates.values()])
    vmin, vmax = float(all_values.min()), float(all_values.max())

    # Individual heatmaps
    print("\nPlotting individual plates...")
    for name, data in plates.items():
        bio, tech = name.split("_")
        title = f"Biological replicate {bio[1]}  ·  Technical replicate {tech}"
        plot_plate_heatmap(data, title, Config.OUTPUT_DIR / f"{name}.png", vmin, vmax)

    # Average heatmap
    print("\nPlotting average plate...")
    average = np.stack(list(plates.values()), axis=0).mean(axis=0)
    plot_plate_heatmap(
        average,
        f"Average plate  (n = {len(plates)})",
        Config.OUTPUT_DIR / "average_plate.png",
        vmin, vmax,
    )

    # Edge-effect analysis
    print("\nEdge-effect correlation analysis...")
    df_all = build_tidy_df(plates)
    plot_edge_correlation(df_all)
    plot_zone_map()
    print_summary(df_all)

    print(f"\nDone – all images saved to {Config.OUTPUT_DIR}")


if __name__ == "__main__":
    main()
