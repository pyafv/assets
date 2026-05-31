import numpy as np
import pyafv as afv


domain_grid = (2, 1)
subdomain_grid = (6, 8)

def main():
    # ========================================================
    # Put MPI setup inside main() so spawned multiprocessing
    # workers do not import/initialize MPI at top level.
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()    # should be 2 for this test
    # ========================================================

    radius = 1.0
    data_multiprocessing = np.load('pts_forces.npz')
    points = data_multiprocessing['points']
    phys = afv.PhysicalParams(r=radius)

    # 2x1 domain decomposition for the two MPI ranks
    if rank == 0:
        domains = afv.decompose_points(points, domain_grid, halo_width=4.01*radius)

    domains = comm.bcast(domains if rank == 0 else None, root=0)

    # Each rank processes its own domain
    domain = domains[rank]
    local_points = domain.local_pts

    # Each rank creates its own simulator with 4 workers
    sim = afv.ParallelFiniteVoronoiSimulator(local_points, phys,
                                             subdomain_grid, n_workers=subdomain_grid[0]*subdomain_grid[1])


    with sim:
        diag = sim.build()

        # gather diagnostics from all ranks
        diag_all = comm.gather(diag, root=0)

        # combine diagnostics on rank 0
        if rank == 0:
            forces = np.zeros_like(points, dtype=float)
            
            for mpi_domain, diag in zip(domains, diag_all):
                owned_local_ids = mpi_domain.owned_local_ids
                owned_global_ids = mpi_domain.local_global_ids[owned_local_ids]

                forces[owned_global_ids] = diag['forces'][owned_local_ids]

            print(f"Forces match or not (MPI)? {np.allclose(forces, data_multiprocessing['diag'])}")  # check forces match
            
            max_diff = np.max(np.abs(forces - data_multiprocessing['diag']))
            print(f'Max difference in forces (MPI): {max_diff}')

    print(f"Rank {rank} finished simulation.")

    if rank == 1:
        sim_serial = afv.FiniteVoronoiSimulator(points, phys)
        diag = sim_serial.build()
        print(f"Forces match or not (serial)? {np.allclose(diag['forces'], data_multiprocessing['diag'])}")  # check forces match
        max_diff_serial = np.max(np.abs(diag['forces'] - data_multiprocessing['diag']))
        print(f'Max difference in forces (serial): {max_diff_serial}')


if __name__ == "__main__":
    import multiprocessing as mp

    mp.set_start_method("spawn", force=True)
    main()
