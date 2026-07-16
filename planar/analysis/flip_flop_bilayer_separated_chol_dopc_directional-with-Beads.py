import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def analyze_flip_flop(resname, headgroup_name, lower_to_upper_filename, upper_to_lower_filename, lower_to_upper_plot, upper_to_lower_plot):
    u = mda.Universe('step7_production.gro', 'concatenated-7us-skipped-10.xtc')
    lipids = u.select_atoms(f'resname {resname} and name {headgroup_name}')

    initial_positions = {}
    for residue in lipids.residues:
        z_position = residue.atoms.select_atoms(f'name {headgroup_name}').positions[0][2]
        initial_positions[residue.resid] = z_position

    times = []
    lower_to_upper_counts = []
    upper_to_lower_counts = []

    for ts in u.trajectory:
        center_of_geometry = lipids.center_of_geometry()
        lower_to_upper_this_frame = 0
        upper_to_lower_this_frame = 0

        for residue in lipids.residues:
            current_z_position = residue.atoms.select_atoms(f'name {headgroup_name}').positions[0][2]
            initial_z_position = initial_positions[residue.resid]

            if current_z_position > center_of_geometry[2] and initial_z_position <= center_of_geometry[2]:
                lower_to_upper_this_frame += 1
                initial_positions[residue.resid] = current_z_position
            elif current_z_position <= center_of_geometry[2] and initial_z_position > center_of_geometry[2]:
                upper_to_lower_this_frame += 1
                initial_positions[residue.resid] = current_z_position

        times.append(ts.time / 1000)
        lower_to_upper_counts.append(lower_to_upper_this_frame)
        upper_to_lower_counts.append(upper_to_lower_this_frame)

    def save_plot_and_data(filename, counts, title, ylabel, plot_filename):
        with open(filename, 'w') as file:
            for time, count in zip(times, counts):
                file.write(f"{time}\t{count}\n")

        plt.figure()
        plt.plot(times, counts)
        plt.xlabel('Time (ns)')
        plt.ylabel(ylabel)
        plt.title(title)
        plt.savefig(plot_filename)
        # plt.show()

    save_plot_and_data(lower_to_upper_filename, lower_to_upper_counts, f'Lower to Upper Leaflet Flip-Flop {resname}', 'Number of Lipids per Frame', lower_to_upper_plot)
    save_plot_and_data(upper_to_lower_filename, upper_to_lower_counts, f'Upper to Lower Leaflet Flip-Flop {resname}', 'Number of Lipids per Frame', upper_to_lower_plot)

analyze_flip_flop('CHOL', 'ROH', 'lower_to_upper_chol.txt', 'upper_to_lower_chol.txt', 'lower_to_upper_chol.png', 'upper_to_lower_chol.png')
analyze_flip_flop('DOPC', 'PO4', 'lower_to_upper_dopc.txt', 'upper_to_lower_dopc.txt', 'lower_to_upper_dopc.png', 'upper_to_lower_dopc.png')

