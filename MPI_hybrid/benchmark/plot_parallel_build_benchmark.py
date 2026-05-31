from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "benchmark_results"
DATA_PATH = RESULTS_DIR / "hpc_parallel_build_times.csv"
OUTPUT_BASENAME = RESULTS_DIR / "hpc_parallel_build_times_nature"
SPEEDUP_DATA_PATH = RESULTS_DIR / "hpc_parallel_build_speedup.csv"

METHOD_ORDER = [
    "serial",
    "parallel_8x6",
    "hybrid_mpi_1x2_parallel_6x4",
    "hybrid_mpi_2x2_parallel_4x3",
    "hybrid_mpi_1x2_parallel_8x6",
    "hybrid_mpi_2x2_parallel_8x6",
]
METHOD_LABELS = {
    "serial": "FiniteVoronoiSimulator",
    "parallel_8x6": "8 x 6",
    "hybrid_mpi_1x2_parallel_6x4": "(1 x 2) x (6 x 4)",
    "hybrid_mpi_2x2_parallel_4x3": "(2 x 2) x (4 x 3)",
    "hybrid_mpi_1x2_parallel_8x6": "(1 x 2) x (8 x 6)",
    "hybrid_mpi_2x2_parallel_8x6": "(2 x 2) x (8 x 6)",
}
COLORS = {
    "serial": "#4A4A4A",
    "parallel_8x6": "#484878",
    "hybrid_mpi_1x2_parallel_6x4": "#F2D6C9",
    "hybrid_mpi_2x2_parallel_4x3": "#E7AA92",
    "hybrid_mpi_1x2_parallel_8x6": "#D77B64",
    "hybrid_mpi_2x2_parallel_8x6": "#B84E45",
}
PANEL_LABEL_EFFECTS = [path_effects.withStroke(linewidth=0.55, foreground="black")]


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Helvetica",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Helvetica",
            "mathtext.cal": "DejaVu Serif Display",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "legend.frameon": False,
        }
    )


def load_data() -> pd.DataFrame:
    paths = [
        RESULTS_DIR / f"{method}_build_times.csv"
        for method in METHOD_ORDER
        if (RESULTS_DIR / f"{method}_build_times.csv").exists()
    ]
    if not paths:
        raise FileNotFoundError(f"no *_build_times.csv files found in {RESULTS_DIR}")

    data = pd.concat((pd.read_csv(path) for path in paths), ignore_index=True)
    data["method"] = pd.Categorical(data["method"], METHOD_ORDER, ordered=True)
    data = data.sort_values(["n_points", "method"])
    data.to_csv(DATA_PATH, index=False)
    return data


def save_pub(fig: mpl.figure.Figure) -> None:
    fig.savefig(f"{OUTPUT_BASENAME}.svg", dpi=300, bbox_inches="tight", transparent=True)
    fig.savefig(f"{OUTPUT_BASENAME}.png", dpi=300, bbox_inches="tight", transparent=True)


def size_label(n_points: int) -> str:
    exponent = int(np.log10(n_points))
    if 10**exponent == n_points:
        return rf"$10^{exponent}$"
    return f"{n_points:,}"


def main() -> None:
    configure_matplotlib()
    data = load_data()
    system_sizes = sorted(data["n_points"].unique())
    methods_present = [method for method in METHOD_ORDER if method in set(data["method"].astype(str))]

    serial = (
        data[data["method"] == "serial"]
        .set_index("n_points")["mean_seconds"]
        .to_dict()
    )
    speedup = data[
        (data["method"] != "serial")
        & data["n_points"].map(lambda n_points: int(n_points) in serial)
    ].copy()
    if not speedup.empty:
        speedup["speedup"] = speedup.apply(
            lambda row: serial[int(row["n_points"])] / row["mean_seconds"],
            axis=1,
        )
    else:
        speedup["speedup"] = np.nan
    speedup.to_csv(SPEEDUP_DATA_PATH, index=False)

    fig = plt.figure(figsize=(7.2, 3.2))
    grid = fig.add_gridspec(1, 2, width_ratios=(1.45, 1.0), wspace=0.24)
    ax_time = fig.add_subplot(grid[0, 0])
    ax_speed = fig.add_subplot(grid[0, 1])

    x = np.arange(len(system_sizes), dtype=float)
    width = 0.115
    offsets = (np.arange(len(METHOD_ORDER)) - 0.5 * (len(METHOD_ORDER) - 1)) * width

    for idx, method in enumerate(METHOD_ORDER):
        method_data = data[data["method"] == method].set_index("n_points")
        values = method_data["mean_seconds"].reindex(system_sizes).to_numpy()
        if np.all(np.isnan(values)):
            continue
        ax_time.bar(
            x + offsets[idx],
            values,
            width=width,
            color=COLORS[method],
            edgecolor="none",
            label=METHOD_LABELS[method],
        )

    ax_time.set_yscale("log")
    ax_time.set_xticks(x)
    ax_time.set_xticklabels([size_label(int(n)) for n in system_sizes])
    ax_time.set_xlabel(r"System size $N$")
    ax_time.set_ylabel("Build time [s]")
    ax_time.grid(axis="y", which="major", color="#D5D5D5", linewidth=0.45)
    ax_time.text(
        -0.12,
        1.04,
        "a",
        transform=ax_time.transAxes,
        fontsize=10,
        fontweight="bold",
        path_effects=PANEL_LABEL_EFFECTS,
        va="bottom",
    )

    for method in METHOD_ORDER[1:]:
        method_data = speedup[speedup["method"] == method].set_index("n_points")
        values = method_data["speedup"].reindex(system_sizes).to_numpy()
        valid = ~np.isnan(values)
        if not np.any(valid):
            continue
        ax_speed.plot(
            np.asarray(system_sizes)[valid],
            values[valid],
            marker="o",
            markersize=3.5,
            linewidth=1.25,
            color=COLORS[method],
            label=METHOD_LABELS[method],
            clip_on=False,
        )

    ax_speed.axhline(1.0, color="#4A4A4A", linewidth=0.7, linestyle="--")
    ax_speed.set_xscale("log")
    ax_speed.set_ylim(bottom=0.0)
    ax_speed.set_xticks(system_sizes)
    ax_speed.set_xticklabels([size_label(int(n)) for n in system_sizes])
    ax_speed.set_xlabel(r"System size $N$")
    ax_speed.set_ylabel("Speedup over FiniteVoronoiSimulator")
    ax_speed.grid(axis="y", color="#D5D5D5", linewidth=0.45)
    ax_speed.text(
        -0.16,
        1.04,
        "b",
        transform=ax_speed.transAxes,
        fontsize=10,
        fontweight="bold",
        path_effects=PANEL_LABEL_EFFECTS,
        va="bottom",
    )

    handles, labels = ax_time.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.52, 1.0),
        ncol=len(methods_present),
        columnspacing=0.9,
        handlelength=1.1,
    )

    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.22, top=0.82, wspace=0.18)
    save_pub(fig)
    print(f"saved {OUTPUT_BASENAME}.svg")
    print(f"saved {DATA_PATH}")
    print(f"saved {SPEEDUP_DATA_PATH}")


if __name__ == "__main__":
    main()
