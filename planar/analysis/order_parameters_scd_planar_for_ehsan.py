import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore', message='Failed to guess the mass for the following atom types:')

u = mda.Universe('step7_production.gro', 'concatenated-7us-skipped-10.xtc')

sn1_atoms = u.select_atoms('name C1A D2A C3A C4A')
sn2_atoms = u.select_atoms('name C1B D2B C3B C4B')

carbon_numbers = list(range(1, 5))

avg_scd_sn1 = []
avg_scd_sn2 = []

for ts in u.trajectory:
    scd_sn1_frame = []
    scd_sn2_frame = []

    for atoms, scd_list in zip([sn1_atoms, sn2_atoms], [scd_sn1_frame, scd_sn2_frame]):
        for i in range(len(atoms) - 1):
            bond_vector = atoms[i+1].position - atoms[i].position
            bond_vector /= np.linalg.norm(bond_vector)
            cos_theta = bond_vector[2]
            scd = 0.5 * (3 * cos_theta**2 - 1)
            scd_list.append(scd)

    avg_scd_sn1.append(np.mean(scd_sn1_frame))
    avg_scd_sn2.append(np.mean(scd_sn2_frame))

final_scd_sn1 = [np.mean([avg_scd_sn1[i] for i in range(j, len(avg_scd_sn1), 4)]) for j in range(4)]
final_scd_sn2 = [np.mean([avg_scd_sn2[i] for i in range(j, len(avg_scd_sn2), 4)]) for j in range(4)]

with open('scd_data.dat', 'w') as file:
    for cn, scd1, scd2 in zip(carbon_numbers, final_scd_sn1, final_scd_sn2):
        file.write(f"{cn}\t{scd1:.3f}\t{scd2:.3f}\n")

plt.plot(carbon_numbers, final_scd_sn1, marker='o', label='Sn1')
plt.plot(carbon_numbers, final_scd_sn2, marker='o', label='Sn2')
plt.xlabel('Carbon Number')
plt.ylabel('Order Parameter S_CD')
plt.title('SCD vs. Carbon Number for Sn1 and Sn2')
plt.legend()
plt.grid(True)
plt.savefig('scd_plot.pdf')
#plt.show()
