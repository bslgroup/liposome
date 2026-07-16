import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def analyze_lipid_distribution(universe_path, trajectory_path, data_filename, plot_filename):
    u = mda.Universe(universe_path, trajectory_path)
    chol = u.select_atoms('resname CHOL')
    dopc = u.select_atoms('resname DOPC')

    def count_lipids_in_leaflets(lipids):
        upper = 0
        lower = 0
        center_of_geometry = lipids.center_of_geometry()[2]
        for lipid in lipids.residues:
            if lipid.atoms.center_of_geometry()[2] > center_of_geometry:
                upper += 1
            else:
                lower += 1
        return upper, lower

    times = []
    chol_upper, chol_lower, dopc_upper, dopc_lower = [], [], [], []
    
    for ts in u.trajectory:
        cu, cl = count_lipids_in_leaflets(chol)
        du, dl = count_lipids_in_leaflets(dopc)
        
        chol_upper.append(cu)
        chol_lower.append(cl)
        dopc_upper.append(du)
        dopc_lower.append(dl)
        
        times.append(ts.time / 1000)  

    chol_ratio = np.array(chol_upper) / np.array(chol_lower, dtype=np.float64)
    dopc_ratio = np.array(dopc_upper) / np.array(dopc_lower, dtype=np.float64)

    with open(data_filename, 'w') as file:
        file.write("Time (ns)\tCHOL Upper\tDOPC Upper\tCHOL Lower\tDOPC Lower\tCHOL Ratio\tDOPC Ratio\n")
        for time, cu, du, cl, dl, cr, dr in zip(times, chol_upper, dopc_upper, chol_lower, dopc_lower, chol_ratio, dopc_ratio):
            file.write(f"{time}\t{cu}\t{du}\t{cl}\t{dl}\t{cr}\t{dr}\n")

    plt.figure()
    plt.plot(times, chol_ratio)
    plt.xlabel('Time (ns)')
    plt.ylabel('CHOL Ratio (Upper/Lower)')
    plt.title('CHOL Ratio vs Time')
    plt.savefig(plot_filename)
    plt.show()

analyze_lipid_distribution('step7_production.gro', 'concatenated-7us-skipped-10.xtc', 'lipid_distribution_data.txt', 'chol_ratio_plot.png')

