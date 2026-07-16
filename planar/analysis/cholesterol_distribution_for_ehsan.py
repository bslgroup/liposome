import MDAnalysis as mda
from MDAnalysis.analysis.density import DensityAnalysis
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings("ignore", message="Failed to guess the mass for the following atom types: ")

sns.set(style='whitegrid')

plt.rcParams.update({
    'font.size': 14,
    'lines.linewidth': 2,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'legend.fontsize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.edgecolor': 'black',
    'axes.linewidth': 1.5,
    'grid.color': 'gray',
    'grid.alpha': 0.5,
    'grid.linestyle': '--'
})

u = mda.Universe('step7_production.gro', 'concatenated-7us-skipped-10.xtc')

cholesterol_sel = 'resname CHOL'
lipid_sel = 'resname DOPC'

def calculate_and_plot_density(universe, axis, output_file, plot_file):
    chol_atoms = universe.select_atoms(cholesterol_sel)
    chol_density = DensityAnalysis(chol_atoms, delta=1.0)
    chol_density.run()
    chol_density_data = chol_density.results['density'].grid

    dopc_atoms = universe.select_atoms(lipid_sel)
    dopc_density = DensityAnalysis(dopc_atoms, delta=1.0)
    dopc_density.run()
    dopc_density_data = dopc_density.results['density'].grid

    axis_map = {'x': 0, 'y': 1, 'z': 2}
    axis_index = axis_map[axis]

    axes_to_sum = [i for i in range(3) if i != axis_index]
    chol_distribution = np.sum(chol_density_data, axis=tuple(axes_to_sum))
    dopc_distribution = np.sum(dopc_density_data, axis=tuple(axes_to_sum))

    chol_distribution /= np.sum(chol_distribution)
    dopc_distribution /= np.sum(dopc_distribution)

    max_length = max(len(chol_distribution), len(dopc_distribution))
    chol_distribution = np.pad(chol_distribution, (0, max_length - len(chol_distribution)), 'constant')
    dopc_distribution = np.pad(dopc_distribution, (0, max_length - len(dopc_distribution)), 'constant')

    bin_edges = np.arange(max_length)

    np.savetxt(output_file, np.column_stack((bin_edges, chol_distribution, dopc_distribution)),
               header=f'Bin edges, CHOL distribution, DOPC distribution along {axis.upper()}-axis')

    plt.figure(figsize=(10, 8))
    plt.plot(bin_edges, chol_distribution, label='CHOL', color='blue')
    plt.plot(bin_edges, dopc_distribution, label='DOPC', color='red')
    plt.xlabel(f'Position along {axis.upper()}-axis (Å)')
    plt.ylabel('Density')  
    plt.title(f'Distribution along {axis.upper()}-axis')
    plt.legend()
    plt.grid(True)
    plt.savefig(plot_file, format='pdf')
    plt.close()

axes = ['x', 'y', 'z']
for axis in axes:
    output_file = f'distribution_{axis}.txt'
    plot_file = f'distribution_{axis}.pdf'
    calculate_and_plot_density(u, axis, output_file, plot_file)

print("Analysis complete. Distribution data and plots saved.")

