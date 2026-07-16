import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('step7_production.gro', 'extract_skipped-100_from_concatenated.xtc')

times = []
outer_radii = []
inner_radii = []
outer_std_devs = []
inner_std_devs = []

for ts in u.trajectory:
    time_ns = u.trajectory.time / 10.0  
    times.append(time_ns)

    sphere_center = u.atoms.center_of_geometry()

    po4_groups = u.select_atoms('resname DOPC and name PO4')
    distances = np.linalg.norm(po4_groups.positions - sphere_center, axis=1)
    median_distance = np.median(distances)
    
    outer_leaflet = po4_groups[distances > median_distance]
    inner_leaflet = po4_groups[distances <= median_distance]

    outer_distances = np.linalg.norm(outer_leaflet.positions - sphere_center, axis=1)
    inner_distances = np.linalg.norm(inner_leaflet.positions - sphere_center, axis=1)

    if ts.frame == 0:
        print(f"Initial frame outer leaflet atoms: {len(outer_leaflet)}")
        print(f"Initial frame inner leaflet atoms: {len(inner_leaflet)}")

    outer_radii.append(np.mean(outer_distances))
    outer_std_devs.append(np.std(outer_distances))

    if len(inner_distances) > 0:
        inner_radii.append(np.mean(inner_distances))
        inner_std_devs.append(np.std(inner_distances))
    else:
        inner_radii.append(np.nan)
        inner_std_devs.append(np.nan)


with open('leaflet_data.txt', 'w') as f:
    f.write('Time (ns)\tOuterRadius\tInnerRadius\tOuterStd Dev\tInnerStdDev\n')
    for i in range(len(times)):
        f.write(f'{times[i]}\t{outer_radii[i]}\t{inner_radii[i]}\t{outer_std_devs[i]}\t{inner_std_devs[i]}\n')

fig, axs = plt.subplots(2, 1, figsize=(10, 8))

axs[0].plot(times, outer_std_devs, label='Outer Leaflet')
axs[0].plot(times, inner_std_devs, label='Inner Leaflet')
axs[0].set_xlabel('Time (ns)')
axs[0].set_ylabel('Standard Deviation')
axs[0].legend()
axs[0].set_title('Smoothness of inner and outer leaflet')

axs[1].plot(times, outer_radii, label='Outer Leaflet')
axs[1].plot(times, inner_radii, label='Inner Leaflet')
axs[1].set_xlabel('Time (ns)')
axs[1].set_ylabel('Radius (Angstrom)')
axs[1].legend()
axs[1].set_title('Radius of inner and outer leaflet')

plt.tight_layout()
plt.savefig('Radius_smoothness_plots.png')
plt.show()

