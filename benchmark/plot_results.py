import numpy as np
import matplotlib.pyplot as plt


def plt_init():
    plt.rcParams.update({
        # requires a working LaTeX installation & slower than Matplotlib's mathtext
        "text.usetex": False,
        "font.family": "Helvetica",    # or other sans-serif font like Arial
        # "font.sans-serif": "DejaVu Sans"
    })
    plt.rcParams["mathtext.fontset"] = "custom"
    # change the math font to Helvetica
    plt.rcParams["mathtext.rm"] = "Helvetica"
    # use the dejavusans version \mathcal{}; set to "DejaVu Sans" for DejaVu normal math font
    plt.rcParams["mathtext.cal"] = "DejaVu Serif Display"

    # save svg with embedded fonts instead of converting text to paths
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["ps.usedistiller"] = "xpdf"    # for scalable eps file

    # depends on whether you want to adjust the distance between subplots by yourself
    plt.rcParams['figure.constrained_layout.use'] = False
    # You can always add bbox_inches='tight' when save figures

    # plt.rcParams['ytick.right'] = True
    # plt.rcParams['xtick.top'] = True
    plt.rc('xtick', labelsize=12, direction='in')
    plt.rc('ytick', labelsize=12, direction='in')
    plt.rc('axes', labelsize=14, titlesize=14)
    plt.rc('legend', fontsize=14)

plt_init()

data_matlab = np.loadtxt('different_N_times_matlab.csv', delimiter=',')
data = np.loadtxt('different_N_times.csv', delimiter=',')
data_Voronoi = np.loadtxt('different_N_times_Voronoi.csv', delimiter=',')

Ns = (10**np.linspace(1, 3, 10)).astype(int)

means_matlab = np.mean(data_matlab, axis=1)
stds_matlab = np.std(data_matlab, axis=1)

means = np.mean(data, axis=1)
stds = np.std(data, axis=1)

means_Voronoi = np.mean(data_Voronoi, axis=1)
stds_Voronoi = np.std(data_Voronoi, axis=1)

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(Ns, means_matlab, marker='o', label='MATLAB', color='C0', clip_on=False)
# ax.fill_between(Ns, means_matlab - stds_matlab, means_matlab + stds_matlab, alpha=0.8)

ax.plot(Ns, means, marker='o', label='PyAFV', color='C3', clip_on=False)
# ax.fill_between(Ns, means - stds, means + stds, alpha=0.8)

ax.plot(Ns, means_Voronoi, marker='o', label='SciPy Voronoi', color='C4', clip_on=False)

xx = np.array([400, 800])
yy = xx/170.
ax.plot(xx, yy, 'k--', lw=1.5)
ax.text(550, 2.05, r'${\sim}N$', fontsize=14)

xx = np.array([300, 600])
yy = xx**1.5/90.
ax.plot(xx, yy, 'k--', lw=1.5)
ax.text(210, 120, r'${\sim}N^{3/2}$', fontsize=14)


ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlabel('Number of cells $N$')
ax.set_ylabel(r'Runtime [$\mathrm{s}$]')
ax.set_title('Benchmark for $10^3$ integration steps')

ax.grid(True, which="both", lw=1, alpha=0.2)

ax.set_xlim(10, 1000)
# ax.set_ylim(bottom=min(means_Voronoi))

plt.legend(fontsize=12.5)

print(Ns)
print(means_matlab/means)

# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)

plt.savefig('runtime_comparison.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()

