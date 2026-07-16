import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def analyze_lipid_distribution_vesicle(universe_path, trajectory_path, data_filename, plot_filename):
    u = mda.Universe(universe_path, trajectory_path)
    chol = u.select_atoms('resname CHOL')
    dopc = u.select_atoms('resname DOPC')
    all_lipids = u.select_atoms('resname CHOL or resname DOPC')

    def median_distance(lipids):
        center_of_vesicle = lipids.center_of_geometry()
        distances = [np.linalg.norm(residue.atoms.center_of_geometry() - center_of_vesicle) for residue in lipids.residues]
        return np.median(distances)

    median_radius = median_distance(all_lipids)

    def count_lipids_in_leaflets(lipids, median_radius, center_of_vesicle):
        upper = 0
        lower = 0
        for lipid in lipids.residues:
            distance = np.linalg.norm(lipid.atoms.center_of_geometry() - center_of_vesicle)
            if distance < median_radius:
                lower += 1
            else:
                upper += 1
        return upper, lower

    times = []
    chol_upper, chol_lower, dopc_upper, dopc_lower = [], [], [], []

    for ts in u.trajectory:
        center_of_vesicle = all_lipids.center_of_geometry()
        cu, cl = count_lipids_in_leaflets(chol, median_radius, center_of_vesicle)
        du, dl = count_lipids_in_leaflets(dopc, median_radius, center_of_vesicle)
        
        chol_upper.append(cu)
        chol_lower.append(cl)
        dopc_upper.append(du)
        dopc_lower.append(dl)
        
        times.append(ts.time / 1000)  

    chol_ratio = np.array(chol_upper, dtype=np.float64) / np.array(chol_lower, dtype=np.float64)
    dopc_ratio = np.array(dopc_upper, dtype=np.float64) / np.array(dopc_lower, dtype=np.float64)

    with open(data_filename, 'w') as file:
        file.write("Time (ns)\tCHOL_Upper\tDOPC_Upper\tCHOL_Lower\tDOPC_Lower\tCHOL_Ratio\tDOPC_Ratio\n")
        for time, cu, du, cl, dl, cr, dr in zip(times, chol_upper, dopc_upper, chol_lower, dopc_lower, chol_ratio, dopc_ratio):
            file.write(f"{time}\t{cu}\t{du}\t{cl}\t{dl}\t{cr}\t{dr}\n")

    plt.figure()
    plt.plot(times, chol_ratio)
    plt.xlabel('Time (ns)')
    plt.ylabel('CHOL Ratio (Upper/Lower)')
    plt.title('CHOL Ratio vs Time in Spherical Vesicle')
    plt.savefig(plot_filename)
    plt.show()

analyze_lipid_distribution_vesicle('step7_production.gro', 'all.xtc', 'lipid_distribution_vesicle_data.txt', 'chol_ratio_vesicle_plot.png')

