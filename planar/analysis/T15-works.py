import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('step7_production.gro', 'all-7micro.xtc')

all_PO4 = u.select_atoms('name PO4')

times = []
thicknesses = []

for ts in u.trajectory[::100]:  # Skipping frames for efficiency
    median_z = np.median(all_PO4.positions[:, 2])

    upper_leaflet = u.select_atoms(f'name PO4 and prop z > {median_z}')
    lower_leaflet = u.select_atoms(f'name PO4 and prop z < {median_z}')

    upper_cog = upper_leaflet.center_of_geometry()
    lower_cog = lower_leaflet.center_of_geometry()

    thickness = abs(upper_cog[2] - lower_cog[2])
    
    times.append(u.trajectory.time / 1000)  # Convert time to nanoseconds
    thicknesses.append(thickness)

with open('membrane_thickness_data.txt', 'w') as file:
    for time, thickness in zip(times, thicknesses):
        file.write(f"{time}\t{thickness}\n")

plt.plot(times, thicknesses)
plt.xlabel('Time (ns)')
plt.ylabel('Membrane Thickness (Å)')
plt.title('Membrane Thickness Over Time')
plt.ylim(0, 45)  # Set y-axis limits
plt.savefig('membrane_thickness_plot.png')
#plt.show()

