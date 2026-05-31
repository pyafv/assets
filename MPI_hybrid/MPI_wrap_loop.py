# MPI_wrap.py

import numpy as np
import pyafv as afv
from tqdm import tqdm
import matplotlib.pyplot as plt


domain_grid = (2, 1)
subdomain_grid = (4, 6)


def main():
    # ========================================================
    # Put MPI setup inside main() so spawned multiprocessing
    # workers do not import/initialize MPI at top level.
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()    # should be 2 for this test
    # ========================================================

    import os

    if hasattr(os, "sched_getaffinity"):
        print(
            f"rank {rank}: affinity size = {len(os.sched_getaffinity(0))}, "
            f"cpus = {sorted(os.sched_getaffinity(0))[:8]}..."
        )

    radius = 1.0
    points = np.random.default_rng(42).random((10_000, 2)) * 100.0
    phys = afv.PhysicalParams(r=radius)

    # 2x1 domain decomposition for the two MPI ranks
    if rank == 0:
        domains = afv.decompose_points(points, domain_grid, halo_width=4.01*radius)

    domains = comm.bcast(domains if rank == 0 else None, root=0)

    # Each rank processes its own domain
    domain = domains[rank]
    local_points = domain.local_pts

    # Each rank creates its own simulator with 4 workers
    sim = afv.ParallelFiniteVoronoiSimulator(local_points, phys, subdomain_grid, n_workers=subdomain_grid[0]*subdomain_grid[1])

    dt = 0.01
    disable_tqdm = (rank != 0)

    with sim:
        for _ in tqdm(range(100), disable=disable_tqdm):
            diag = sim.build()

            # gather diagnostics from all ranks
            diag_all = comm.gather(diag, root=0)

            # combine diagnostics on rank 0
            if rank == 0:
                forces = np.zeros_like(points, dtype=float)
                areas = np.zeros(points.shape[0], dtype=float)
                perimeters = np.zeros(points.shape[0], dtype=float)
                # ... add more diagnostics as needed ...

                for mpi_domain, diag in zip(domains, diag_all):
                    owned_local_ids = mpi_domain.owned_local_ids
                    owned_global_ids = mpi_domain.local_global_ids[owned_local_ids]

                    forces[owned_global_ids] = diag["forces"][owned_local_ids]
                    areas[owned_global_ids] = diag["areas"][owned_local_ids]
                    perimeters[owned_global_ids] = diag["perimeters"][owned_local_ids]

                diag_combined = {
                    "forces": forces,
                    "areas": areas,
                    "perimeters": perimeters,
                }
                

                # Update points
                points += forces * dt

                # upate domain decomposition with new points
                domains = afv.decompose_points(points, domain_grid, halo_width=4.01*radius)

            domains = comm.bcast(domains, root=0)
            
            domain = domains[rank]

            local_points = domain.local_pts
            sim.update_positions(local_points)
        
        diag = sim.build(plot_mode=True)
        diag_all = comm.gather(diag, root=0)

        if rank == 0:      # use only one rank to plot
            fig, ax = plt.subplots()
            for idx in range(size):
                diag = diag_all[idx]
                domain = domains[idx]
                local_points = domain.local_pts

                afv.visualize_2d_parallel(local_points, diag, r=radius, ax=ax,
                        selected=domain.owned_local_ids)  # remove halo points for each rank's domain

            ax.set_xlim(-10, 110)
            ax.set_ylim(-10, 110)
            plt.savefig("mpi_parallel_voronoi.png", dpi=300)

    print(f"\nRank {rank} finished simulation.")


if __name__ == "__main__":
    import multiprocessing as mp

    mp.set_start_method("spawn", force=True)
    main()

