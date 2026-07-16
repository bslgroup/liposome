import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def analyze_flip_flop(resname, data_filename, plot_filename):
    u = mda.Universe('step7_production.gro', 'concatenated-7us-skipped-10.xtc')
    lipids = u.select_atoms(f'resname {resname}')  

    initial_positions = {}
    for residue in lipids.residues:
        center_of_geometry = lipids.center_of_geometry()
        z_distance = residue.atoms.center_of_geometry()[2] - center_of_geometry[2]
        initial_positions[residue.resid] = z_distance

    times = []
    flip_flop_counts = []

    for ts in u.trajectory:
        center_of_geometry = lipids.center_of_geometry()
        flip_flops_this_frame = 0
        for residue in lipids.residues:
            current_z_distance = residue.atoms.center_of_geometry()[2] - center_of_geometry[2]
            if np.sign(current_z_distance) != np.sign(initial_positions[residue.resid]):
                flip_flops_this_frame += 1
                initial_positions[residue.resid] = current_z_distance

        times.append(ts.time / 1000)  
        flip_flop_counts.append(flip_flops_this_frame)

    with open(data_filename, 'w') as file:
        for time, count in zip(times, flip_flop_counts):
            file.write(f"{time}\t{count}\n")

    plt.figure()  
    plt.plot(times, flip_flop_counts)
    plt.xlabel('Time (ns)')
    plt.ylabel('Number of Flip-Flopping Residues per Frame')
    plt.title(f'Number of Flip-Flopping {resname} Residues per Time Frame')
    plt.savefig(plot_filename)
    #plt.show()

analyze_flip_flop('CHOL', 'flip_flop_data_chol.txt', 'flip_flop_plot_chol.png')

analyze_flip_flop('DOPC', 'flip_flop_data_dopc.txt', 'flip_flop_plot_dopc.png')

