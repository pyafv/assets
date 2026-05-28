import numpy as np
import tqdm
import pyafv
import matplotlib.pyplot as plt

points = np.random.default_rng(42).random((10_000, 2)) * 100.0
phys = pyafv.PhysicalParams(r=1.0)

sim = pyafv.ParallelFiniteVoronoiSimulator(
    points,
    phys,
    grid_shape=(2, 2),
    n_workers=4,
)

diag = sim.build(plot_mode=True)

print(diag["pids"])

pyafv.visualize_2d_parallel(points, diag, r=1.0)

plt.savefig('test.png', dpi=300, bbox_inches='tight')
