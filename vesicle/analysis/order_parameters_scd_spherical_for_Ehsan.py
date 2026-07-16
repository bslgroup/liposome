# Author: Mortaza Derakhshani-Molayousefi 	University of Arkansas		Email: derakhshanimortaza@gmail.com
# This script computes and plots the order parameter (SCD) for Spherical Lipid bilayer 

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore") # , message="Warning Message")
def calculate_order_parameter(atom1, atom2, bilayer_center):
    bond_vector = atom2.position - atom1.position
    bond_vector /= np.linalg.norm(bond_vector)
    normal_vector = atom1.position - bilayer_center
    normal_vector /= np.linalg.norm(normal_vector)
    cos_theta = np.dot(bond_vector, normal_vector)
    scd = 0.5 * (3 * cos_theta**2 - 1)
    return scd

u = mda.Universe('step7_production.gro', 'extract_skipped-100_from_concatenated.xtc')
bilayer_center = u.select_atoms('name C1A C1B').center_of_mass()

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
            scd = calculate_order_parameter(atoms[i], atoms[i+1], bilayer_center)
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
plt.show()

