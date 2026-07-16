import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

u = mda.Universe('step7_production.gro', 'concatenated-7us-skipped-10.xtc')

cholesterols = u.select_atoms('name ROH')

num_cholesterols_per_leaflet = len(cholesterols) / 2  

times = []
areas_per_frame = []

for ts in u.trajectory:
    box_area = ts.dimensions[0] * ts.dimensions[1]  
    area_per_cholesterol = box_area / num_cholesterols_per_leaflet
    times.append(ts.time / 1000)  
    areas_per_frame.append(area_per_cholesterol)

with open('cholesterol_areas.txt', 'w') as f:
    for time, area in zip(times, areas_per_frame):
        f.write(f"{time}\t{area}\n")

plt.plot(times, areas_per_frame)
plt.xlabel('Time (ns)')
plt.ylabel('Average Area per Cholesterol (Å²)')
plt.title('Average Cholesterol Area')
plt.savefig('chol_area_plot.png')
#plt.show()

