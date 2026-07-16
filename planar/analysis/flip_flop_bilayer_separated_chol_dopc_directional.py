import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def analyze_flip_flop(resname):
    u = mda.Universe('step7_production.gro', 'all.xtc')
    lipids = u.select_atoms(f'resname {resname}')

    initial_positions = {}
    for residue in lipids.residues:
        center_of_geometry = lipids.center_of_geometry()
        z_distance = residue.atoms.center_of_geometry()[2] - center_of_geometry[2]
        initial_positions[residue.resid] = z_distance

    lower_to_upper = []
    upper_to_lower = []
    times = []

    for ts in u.trajectory:
        center_of_geometry = lipids.center_of_geometry()
        lower_to_upper_this_frame = 0
        upper_to_lower_this_frame = 0
        for residue in lipids.residues:
            current_z_distance = residue.atoms.center_of_geometry()[2] - center_of_geometry[2]
            initial_z_distance = initial_positions[residue.resid]

            if np.sign(current_z_distance) != np.sign(initial_z_distance):
                if np.sign(initial_z_distance) < 0:  
                    lower_to_upper_this_frame += 1
                else:  
                    upper_to_lower_this_frame += 1
                initial_positions[residue.resid] = current_z_distance

        times.append(ts.time / 1000)
        lower_to_upper.append(lower_to_upper_this_frame)
        upper_to_lower.append(upper_to_lower_this_frame)

    data_filenames = [f'{resname}_lower_to_upper.txt', f'{resname}_upper_to_lower.txt']
    plot_filenames = [f'{resname}_lower_to_upper_plot.png', f'{resname}_upper_to_lower_plot.png']
    counts = [lower_to_upper, upper_to_lower]
    directions = ['Lower to Upper', 'Upper to Lower']

    for data_filename, plot_filename, count, direction in zip(data_filenames, plot_filenames, counts, directions):
        with open(data_filename, 'w') as file:
            for time, c in zip(times, count):
                file.write(f"{time}\t{c}\n")

        plt.figure()
        plt.plot(times, count)
        plt.xlabel('Time (ns)')
        plt.ylabel(f'Number of {direction} Flip-Flops')
        plt.title(f'Number of {direction} {resname} Flip-Flops vs. Time (ns)')
        plt.savefig(plot_filename)

analyze_flip_flop('CHOL')
analyze_flip_flop('DOPC')

