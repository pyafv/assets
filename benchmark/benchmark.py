import numpy as np
import pyafv as afv
import time

pts = np.loadtxt('init_pts.csv', delimiter=',')
phys = afv.PhysicalParams()

sim = afv.FiniteVoronoiSimulator(pts, phys)

dt = 0.01
T = 30.
num_steps = int(T / dt)

start = time.time()
for _ in range(num_steps):
    forces = sim.build()['forces']
    pts += forces * dt
    sim.update_positions(pts)

end = time.time()
print(f"Simulation completed in {end - start:.6f} seconds.")

sim.plot_2d(show=True)
