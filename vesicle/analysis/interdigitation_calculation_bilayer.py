import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('step7_production.gro', 'all.xtc')

all_tails = u.select_atoms('name GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B')  

times = []
interdigitation_values = []

for ts in u.trajectory[::1]:  
    median_z = np.median(all_tails.positions[:, 2])

    upper_tails = u.select_atoms(f'name GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B and prop z > {median_z}')
    lower_tails = u.select_atoms(f'name GL1 GL2 C1A D2A C3A C4A C1B D2B C3B C4B and prop z < {median_z}')

    upper_range = upper_tails.positions[:, 2].max() - upper_tails.positions[:, 2].min()
    lower_range = lower_tails.positions[:, 2].max() - lower_tails.positions[:, 2].min()

    interdigitation = abs(upper_range - lower_range)
    
    times.append(ts.time / 1000)
    interdigitation_values.append(interdigitation)

with open('interdigitation_data.txt', 'w') as file:
    for time, interdigitation in zip(times, interdigitation_values):
        file.write(f"{time}\t{interdigitation}\n")

plt.plot(times, interdigitation_values)
plt.xlabel('Time (ns)')
plt.ylabel('Interdigitation (Å)')
plt.title('Interdigitation Over Time')
plt.savefig('interdigitation_plot.png')
plt.show()

