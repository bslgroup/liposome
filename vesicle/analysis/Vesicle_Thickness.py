import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('step7_production.gro', 'extract_skipped-100_from_concatenated.xtc')
all_PO4 = u.select_atoms('name PO4')

times = []
thicknesses = []

for ts in u.trajectory[::1]:  # Skipping frames for efficiency
    vesicle_center = all_PO4.center_of_geometry()

    radial_distances = np.linalg.norm(all_PO4.positions - vesicle_center, axis=1)

    median_distance = np.median(radial_distances)
    inner_leaflet = all_PO4[radial_distances < median_distance]
    outer_leaflet = all_PO4[radial_distances >= median_distance]

    inner_distances = np.linalg.norm(inner_leaflet.positions - vesicle_center, axis=1)
    outer_distances = np.linalg.norm(outer_leaflet.positions - vesicle_center, axis=1)
    avg_thickness = np.mean(outer_distances) - np.mean(inner_distances)

    times.append(ts.time / 10)  # Convert time to nanoseconds
    thicknesses.append(avg_thickness)

with open('vesicle_thickness_data.txt', 'w') as file:
    for time, thickness in zip(times, thicknesses):
        file.write(f"{time}\t{thickness}\n")

plt.plot(times, thicknesses)
plt.xlabel('Time (ns)')
plt.ylabel('Vesicle Thickness (Å)')
plt.title('Vesicle Thickness Over Time')
plt.savefig('vesicle_thickness_plot.png')
plt.show()

~                
