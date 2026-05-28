from __future__ import annotations

import csv
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pyafv as afv


METHOD = "parallel_2x1"
GRID_SHAPE = (2, 1)
N_WORKERS = 2
SYSTEM_SIZES = [100, 1_000, 10_000, 100_000, 1_000_000]
REPEATS = 10
WARMUP_BUILDS = 3
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


def time_builds(sim, point_sets: list[np.ndarray]) -> tuple[float, float]:
    sim.update_positions(point_sets[0])
    for _ in range(WARMUP_BUILDS):
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
    csv_path = OUTPUT_DIR / f"{METHOD}_build_times.csv"
    grid_label = f"{GRID_SHAPE[0]}x{GRID_SHAPE[1]}"

    phys = afv.PhysicalParams(r=RADIUS)
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["method", "grid_shape", "n_workers", "n_points", "mean_seconds", "std_seconds"])

        for size_idx, n_points in enumerate(SYSTEM_SIZES):
            print(f"{METHOD}: N={n_points:,}", flush=True)
            point_sets = make_point_sets(n_points, size_idx)
            sim = afv.ParallelFiniteVoronoiSimulator(
                point_sets[0],
                phys,
                grid_shape=GRID_SHAPE,
                n_workers=N_WORKERS,
            )
            with sim:
                mean_time, std_time = time_builds(sim, point_sets)
            writer.writerow([METHOD, grid_label, N_WORKERS, n_points, mean_time, std_time])
            handle.flush()
            print(f"{METHOD}: N={n_points:,} mean={mean_time:.6g}s std={std_time:.6g}s", flush=True)

    print(f"saved {csv_path}", flush=True)


if __name__ == "__main__":
    main()
