import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

u = mda.Universe('step7_production.gro', 'extract_skipped-100_from_concatenated.xtc')

cholesterols = u.select_atoms('resname CHOL')  

num_cholesterols = len(cholesterols.residues) // 2  

areas_per_cholesterol = []
radii = []
times = []

for ts in u.trajectory:
    center_of_geometry = cholesterols.center_of_geometry()

    distances = np.linalg.norm(cholesterols.positions - center_of_geometry, axis=1)

    vesicle_radius = np.mean(distances)
    radii.append(vesicle_radius)

    vesicle_surface_area = 4 * np.pi * vesicle_radius**2

    area_per_cholesterol = vesicle_surface_area / num_cholesterols
    areas_per_cholesterol.append(area_per_cholesterol)

    times.append(ts.time / 10)

average_radius = np.mean(radii)

with open('area_per_cholesterol_data.txt', 'w') as file:
    for time, area in zip(times, areas_per_cholesterol):
        file.write(f"{time}\t{area}\n")

plt.plot(times, areas_per_cholesterol)
plt.xlabel('Time (ns)')
plt.ylabel('Average Area per Cholesterol (Å²)')
plt.title('Average Area per Cholesterol')
plt.savefig('area_per_cholesterol_plot.png')
#plt.show()


