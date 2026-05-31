from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pyafv as afv


METHOD = "hybrid_mpi_1x2_parallel_6x4"
DOMAIN_GRID = (1, 2)
SUBDOMAIN_GRID = (6, 4)
N_WORKERS = SUBDOMAIN_GRID[0] * SUBDOMAIN_GRID[1]
SYSTEM_SIZES = [100, 1_000, 10_000, 100_000, 1_000_000]
REPEATS = 10
WARMUP_BUILDS = 3
SEED = 42
RADIUS = 1.0
HALO_WIDTH = 4.01 * RADIUS
OUTPUT_DIR = Path(__file__).resolve().parent / "benchmark_results"


def make_points(n_points: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    box_size = np.sqrt(n_points) * RADIUS
    return rng.random((n_points, 2)) * box_size


def make_initial_local_points(comm, rank: int, n_points: int, seed: int) -> np.ndarray:
    if rank == 0:
        points = make_points(n_points, seed)
        domains = afv.decompose_points(points, DOMAIN_GRID, halo_width=HALO_WIDTH)
    else:
        domains = None

    domains = comm.bcast(domains, root=0)
    return domains[rank].local_pts


def run_step(comm, rank: int, sim, points: np.ndarray | None, n_points: int) -> None:
    if rank == 0:
        domains = afv.decompose_points(points, DOMAIN_GRID, halo_width=HALO_WIDTH)
    else:
        domains = None

    domains = comm.bcast(domains, root=0)
    domain = domains[rank]

    sim.update_positions(domain.local_pts)
    diag = sim.build(connect=False)

    owned_local_ids = domain.owned_local_ids
    owned_global_ids = domain.local_global_ids[owned_local_ids]
    local_forces = diag["forces"][owned_local_ids]

    gathered = comm.gather((owned_global_ids, local_forces), root=0)

    if rank == 0:
        forces = np.zeros((n_points, 2), dtype=float)
        for global_ids, rank_forces in gathered:
            forces[global_ids] = rank_forces


def time_steps(comm, MPI, rank: int, sim, n_points: int, size_idx: int) -> tuple[float | None, float | None]:
    for warmup_idx in range(WARMUP_BUILDS):
        seed = SEED + 1_000_000 + size_idx * WARMUP_BUILDS + warmup_idx
        points = make_points(n_points, seed) if rank == 0 else None
        run_step(comm, rank, sim, points, n_points)

    times = np.empty(REPEATS, dtype=float) if rank == 0 else None
    for repeat_idx in range(REPEATS):
        seed = SEED + size_idx * REPEATS + repeat_idx
        points = make_points(n_points, seed) if rank == 0 else None

        comm.Barrier()
        t0 = MPI.Wtime()
        run_step(comm, rank, sim, points, n_points)
        local_time = MPI.Wtime() - t0

        global_time = comm.reduce(local_time, op=MPI.MAX, root=0)
        if rank == 0:
            times[repeat_idx] = global_time

    if rank == 0:
        return float(np.mean(times)), float(np.std(times))
    return None, None


def main() -> None:
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    n_domains = DOMAIN_GRID[0] * DOMAIN_GRID[1]
    if size != n_domains:
        raise RuntimeError(f"run with {n_domains} MPI ranks for DOMAIN_GRID={DOMAIN_GRID}")

    if hasattr(os, "sched_getaffinity"):
        cpus = sorted(os.sched_getaffinity(0))
        print(f"rank {rank}: affinity size={len(cpus)}, cpus={cpus[:8]}...", flush=True)

    phys = afv.PhysicalParams(r=RADIUS)
    csv_path = OUTPUT_DIR / f"{METHOD}_build_times.csv"

    if rank == 0:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        handle = csv_path.open("w", newline="")
        writer = csv.writer(handle)
        writer.writerow(
            [
                "method",
                "domain_grid",
                "subdomain_grid",
                "n_workers",
                "n_mpi_ranks",
                "n_points",
                "mean_seconds",
                "std_seconds",
            ]
        )
    else:
        handle = None
        writer = None

    for size_idx, n_points in enumerate(SYSTEM_SIZES):
        if rank == 0:
            print(f"{METHOD}: N={n_points:,}", flush=True)

        local_points0 = make_initial_local_points(
            comm,
            rank,
            n_points,
            SEED + 10_000_000 + size_idx,
        )
        sim = afv.ParallelFiniteVoronoiSimulator(
            local_points0,
            phys,
            grid_shape=SUBDOMAIN_GRID,
            n_workers=N_WORKERS,
        )

        with sim:
            mean_time, std_time = time_steps(comm, MPI, rank, sim, n_points, size_idx)

        if rank == 0:
            writer.writerow(
                [
                    METHOD,
                    f"{DOMAIN_GRID[0]}x{DOMAIN_GRID[1]}",
                    f"{SUBDOMAIN_GRID[0]}x{SUBDOMAIN_GRID[1]}",
                    N_WORKERS,
                    size,
                    n_points,
                    mean_time,
                    std_time,
                ]
            )
            handle.flush()
            print(f"{METHOD}: N={n_points:,} mean={mean_time:.6g}s std={std_time:.6g}s", flush=True)

    if rank == 0:
        handle.close()
        print(f"saved {csv_path}", flush=True)


if __name__ == "__main__":
    import multiprocessing as mp

    mp.set_start_method("spawn", force=True)
    main()
