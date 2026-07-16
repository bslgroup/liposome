import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('step7_production.gro', 'concatenated-7us-skipped-10.xtc')

lipids = u.select_atoms('resname DOPC CHOL')  

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

with open('flip_flop_data.txt', 'w') as file:
    for time, count in zip(times, flip_flop_counts):
        file.write(f"{time}\t{count}\n")

plt.plot(times, flip_flop_counts)
plt.xlabel('Time (ns)')
plt.ylabel('Number of Flip-Flopping Residues per Frame')
plt.title('Number of Flip-Flopping Residues per Time Frame')
plt.savefig('flip_flop_plot.png')
plt.show()

