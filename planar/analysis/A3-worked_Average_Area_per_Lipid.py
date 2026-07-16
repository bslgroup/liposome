import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('step7_production.gro', 'concatenated-7us-skipped-10.xtc')

lipids = u.select_atoms('name PO4')

num_lipids_per_leaflet = len(lipids) / 2  

times = []
areas_per_frame = []

for ts in u.trajectory:
    box_area = ts.dimensions[0] * ts.dimensions[1]  # x-size * y-size of the box
    area_per_lipid = box_area / num_lipids_per_leaflet
    times.append(ts.time / 1000)  # Convert time to ns
    areas_per_frame.append(area_per_lipid)

with open('lipid_areas.txt', 'w') as f:
    for time, area in zip(times, areas_per_frame):
        f.write(f"{time}\t{area}\n")

plt.plot(times, areas_per_frame)
plt.xlabel('Time (ns)')
plt.ylabel('Average Area per Lipid (Å²)')
plt.title('Average Lipid Area Over Time')
plt.ylim(60, 80)  # Set y-axis limits
plt.savefig('lipid_area_plot.png')
plt.show()

