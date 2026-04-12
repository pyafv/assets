import numpy as np
import pyafv as afv

phys = afv.PhysicalParams(delta=0.0)

# test of initial forces
pts = np.loadtxt('init_pts.csv', delimiter=',')
sim = afv.FiniteVoronoiSimulator(pts, phys)
forces = sim.build()['forces']

forces_ref = np.loadtxt('init_forces.csv', delimiter=',')
assert np.allclose(forces, forces_ref, atol=1e-8), "Initial forces do not match reference."
print("Initial forces match the reference data within tolerance.")

# test of final forces after simulation
pts = np.loadtxt('final_pts.csv', delimiter=',')
sim.update_positions(pts)
forces = sim.build()['forces']
forces_ref = np.loadtxt('final_forces.csv', delimiter=',')
assert np.allclose(forces, forces_ref, atol=1e-7), "Final forces do not match reference."
print("Final forces match the reference data within tolerance.")

sim.plot_2d(show=True)
