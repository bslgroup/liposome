import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def analyze_flip_flop_vesicle(resname, headgroup_name, prefix):
    u = mda.Universe('step7_production.gro', 'extract_skipped-100_from_concatenated.xtc')
    lipids = u.select_atoms(f'resname {resname} and name {headgroup_name}')

    initial_positions = {}
    center_of_vesicle = lipids.center_of_geometry()
    for residue in lipids.residues:
        headgroup_position = residue.atoms.select_atoms(f'name {headgroup_name}').positions[0]
        initial_positions[residue.resid] = np.linalg.norm(headgroup_position - center_of_vesicle)

    median_radius = np.median(list(initial_positions.values()))

    times = []
    inner_to_outer_counts = []
    outer_to_inner_counts = []

    for ts in u.trajectory:
        center_of_vesicle = lipids.center_of_geometry()
        inner_to_outer_this_frame = 0
        outer_to_inner_this_frame = 0
        for residue in lipids.residues:
            headgroup_position = residue.atoms.select_atoms(f'name {headgroup_name}').positions[0]
            current_radius = np.linalg.norm(headgroup_position - center_of_vesicle)
            if initial_positions[residue.resid] < median_radius and current_radius > median_radius:
                inner_to_outer_this_frame += 1
                initial_positions[residue.resid] = current_radius
            elif initial_positions[residue.resid] > median_radius and current_radius < median_radius:
                outer_to_inner_this_frame += 1
                initial_positions[residue.resid] = current_radius

        times.append(ts.time / 10)
        inner_to_outer_counts.append(inner_to_outer_this_frame)
        outer_to_inner_counts.append(outer_to_inner_this_frame)

    data_filenames = [f'{prefix}_inner_to_outer.txt', f'{prefix}_outer_to_inner.txt']
    plot_filenames = [f'{prefix}_inner_to_outer_plot.png', f'{prefix}_outer_to_inner_plot.png']
    count_data = [inner_to_outer_counts, outer_to_inner_counts]

    for data_filename, plot_filename, counts in zip(data_filenames, plot_filenames, count_data):
        with open(data_filename, 'w') as file:
            for time, count in zip(times, counts):
                file.write(f"{time}\t{count}\n")

        plt.figure()
        plt.plot(times, counts)
        plt.xlabel('Time (ns)')
        plt.ylabel('Number of Flip-Flopping Lipids')
        plt.title(f'{plot_filename[:-4].replace("_", " ").capitalize()}')
        plt.savefig(plot_filename)
        # plt.show()

analyze_flip_flop_vesicle('CHOL', 'ROH', 'chol')
analyze_flip_flop_vesicle('DOPC', 'PO4', 'dopc')

