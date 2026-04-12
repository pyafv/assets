import numpy as np
import time
from scipy.spatial import Voronoi

np.random.seed(42)

Ns = (10**np.linspace(1, 3, 10)).astype(int)

all_pts = []
for N in Ns:
    pts = np.random.rand(N, 2)
    all_pts.append(pts)


all_times = []
for N in Ns:
    print(f"Running simulation for N={N}")
    pts = all_pts[np.where(Ns == N)[0][0]]

    dt = 0.01
    T = 10.
    num_steps = int(T / dt)

    times_N = []
    for _ in range(3):  # Run 3 times for each N to average time
        start = time.time()
        for _ in range(num_steps):
            vor = Voronoi(pts)
            pts = np.random.rand(N, 2)
            # pts += forces * dt

        end = time.time()
        time_taken = end - start
        times_N.append(time_taken)

    all_times.append(times_N)

# Save pts and times separately
# np.savetxt('different_N_pts.csv', np.array(Ns).reshape(-1, 1), delimiter=',')
np.savetxt('different_N_times_Voronoi.csv', all_times, delimiter=',')
