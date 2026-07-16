import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def analyze_flip_flop_vesicle(resname, data_filename, plot_filename):
    u = mda.Universe('step7_production.gro', 'extract_skipped-100_from_concatenated.xtc')
    lipids = u.select_atoms(f'resname {resname}')

    initial_positions = {}
    center_of_vesicle = lipids.center_of_geometry()
    for residue in lipids.residues:
        initial_positions[residue.resid] = np.linalg.norm(residue.atoms.center_of_geometry() - center_of_vesicle)

    median_radius = np.median(list(initial_positions.values()))

    times = []
    flip_flop_counts = []

    for ts in u.trajectory:
        center_of_vesicle = lipids.center_of_geometry()
        flip_flops_this_frame = 0
        for residue in lipids.residues:
            current_radius = np.linalg.norm(residue.atoms.center_of_geometry() - center_of_vesicle)
            if (initial_positions[residue.resid] < median_radius and current_radius > median_radius) or \
               (initial_positions[residue.resid] > median_radius and current_radius < median_radius):
                flip_flops_this_frame += 1
                initial_positions[residue.resid] = current_radius

        times.append(ts.time / 10)
        flip_flop_counts.append(flip_flops_this_frame)

    with open(data_filename, 'w') as file:
        for time, count in zip(times, flip_flop_counts):
            file.write(f"{time}\t{count}\n")

    plt.figure()  
    plt.plot(times, flip_flop_counts)
    plt.xlabel('Time (ns)')
    plt.ylabel('Number of Flip-Flopping Residues per Frame')
    plt.title(f'Number of Flip-Flopping {resname} Residues per Time Frame in Spherical Vesicle')
    plt.savefig(plot_filename)
#    plt.show()

analyze_flip_flop_vesicle('CHOL', 'flip_flop_data_chol.txt', 'flip_flop_plot_chol.png')

analyze_flip_flop_vesicle('DOPC', 'flip_flop_data_dopc.txt', 'flip_flop_plot_dopc.png')

