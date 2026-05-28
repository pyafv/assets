from __future__ import annotations

import csv
import sys
from pathlib import Path
from time import perf_counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pyafv as afv


SYSTEM_SIZES = [100, 1_000, 10_000, 100_000, 1_000_000]
PARALLEL_SETUPS = [
    ("parallel_2x1", (2, 1), 2),
    ("parallel_2x2", (2, 2), 4),
    ("parallel_3x2", (3, 2), 6),
    ("parallel_3x3", (3, 3), 9),
    ("parallel_4x3", (4, 3), 12),
]
REPEATS = 10
SEED = 42
RADIUS = 1.0
OUTPUT_DIR = Path(__file__).resolve().parent / "benchmark_results"


def make_points(n_points: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    box_size = np.sqrt(n_points) * RADIUS
    return rng.random((n_points, 2)) * box_size


def make_point_sets(n_points: int, size_idx: int) -> list[np.ndarray]:
    return [
        make_points(n_points, SEED + size_idx * REPEATS + repeat_idx)
        for repeat_idx in range(REPEATS)
    ]


def time_builds(sim, point_sets: list[np.ndarray], warmup: bool = False) -> tuple[float, float]:
    if warmup:
        sim.update_positions(point_sets[0])
        for _ in range(3):
            sim.build(connect=False)

    times = np.empty(len(point_sets), dtype=float)
    for idx, points in enumerate(point_sets):
        sim.update_positions(points)
        t0 = perf_counter()
        sim.build(connect=False)
        times[idx] = perf_counter() - t0
    return float(np.mean(times)), float(np.std(times))


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIR / "parallel_build_times.csv"
    plot_path = OUTPUT_DIR / "parallel_build_times.png"

    phys = afv.PhysicalParams(r=RADIUS)
    rows = []

    for size_idx, n_points in enumerate(SYSTEM_SIZES):
        point_sets = make_point_sets(n_points, size_idx)

        print(f"N={n_points:,}: serial")
        sim = afv.FiniteVoronoiSimulator(point_sets[0], phys)
        mean_time, std_time = time_builds(sim, point_sets)
        rows.append(["serial", "", 1, n_points, mean_time, std_time])

        for name, grid_shape, n_workers in PARALLEL_SETUPS:
            print(f"N={n_points:,}: {name}")
            sim = afv.ParallelFiniteVoronoiSimulator(
                point_sets[0],
                phys,
                grid_shape=grid_shape,
                n_workers=n_workers,
            )
            with sim:
                mean_time, std_time = time_builds(sim, point_sets, warmup=True)
            rows.append([name, f"{grid_shape[0]}x{grid_shape[1]}", n_workers, n_points, mean_time, std_time])

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "grid_shape", "n_workers", "n_points", "mean_seconds", "std_seconds"])
        writer.writerows(rows)

    methods = ["serial"] + [name for name, _, _ in PARALLEL_SETUPS]
    labels = ["serial", "2x1=2", "2x2=4", "3x2=6", "3x3=9", "4x3=12"]
    x = np.arange(len(SYSTEM_SIZES), dtype=float)
    width = 0.13

    fig, ax = plt.subplots(figsize=(9, 5))
    for idx, (method, label) in enumerate(zip(methods, labels)):
        values = [
            row[4]
            for n_points in SYSTEM_SIZES
            for row in rows
            if row[0] == method and row[3] == n_points
        ]
        offsets = x + (idx - 0.5 * (len(methods) - 1)) * width
        ax.bar(offsets, values, width=width, label=label)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{n_points:,}" for n_points in SYSTEM_SIZES])
    ax.set_xlabel("System size N")
    ax.set_ylabel("Mean build time (s)")
    ax.legend(title="setup")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=300)

    print(f"saved {csv_path}")
    print(f"saved {plot_path}")


if __name__ == "__main__":
    main()
