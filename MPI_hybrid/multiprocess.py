import numpy as np
import pyafv as afv
from tqdm import tqdm
import matplotlib.pyplot as plt


def main():
    radius = 1.0
    points = np.random.default_rng(42).random((100_000, 2)) * np.sqrt(100_000)
    phys = afv.PhysicalParams(r=radius)

    sim = afv.ParallelFiniteVoronoiSimulator(points, phys, (6, 5), n_workers=30)

    dt = 0.01
    with sim:
        for _ in tqdm(range(1000)):
            diag = sim.build()
            points += diag['forces'] * dt
            sim.update_positions(points)

        diag = sim.build(plot_mode=True)

    np.savez('pts_forces.npz', points=points, diag=diag['forces'])

    fig, ax = plt.subplots()
    afv.visualize_2d_parallel(points, diag, r=radius, ax=ax)
    plt.savefig('N_10_5.png', dpi=300)



if __name__ == '__main__':
    main()