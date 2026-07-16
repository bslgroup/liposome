#!/usr/bin/env python

import sys
import os
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import Voronoi
from shapely.geometry import Polygon
from collections import defaultdict

FOLDERS = ["0", "10", "20", "30", "40"]
GRO_FILE = "step7_production.gro"
XTC_FILE = "step7_production.xtc"

CHOL_RESNAME = "CHOL"
CHOL_HEAD_ATOM = "ROH"
DOPC_RESNAME = "DOPC"
DOPC_HEAD_ATOM = "PO4"

NM2_to_ANG2 = 100.0

PX = 1
PY = 1

FRAME_STRIDE = 1

def replicate_points_2d(pts, box_x, box_y, px=1, py=1):
    replicated = []
    for i in range(-px, px+1):
        for j in range(-py, py+1):
            shift = np.array([i*box_x, j*box_y])
            replicated.append(pts + shift)
    return np.concatenate(replicated, axis=0)

def voronoi_2d_pbc(pts, box_x, box_y, px=1, py=1):
    big_pts = replicate_points_2d(pts, box_x, box_y, px=px, py=py)
    vor = Voronoi(big_pts, qhull_options="Qbb Qx")
    return vor

def clip_voronoi_polygons_to_box(vor, pts, box_x, box_y, px=1, py=1):
    from shapely.ops import unary_union

    clip_box = Polygon([(0,0), (box_x,0), (box_x,box_y), (0,box_y)])
    
    coord_to_region = defaultdict(list)
    for i, p in enumerate(vor.points):
        key = (round(p[0],4), round(p[1],4))
        coord_to_region[key].append(vor.point_region[i])

    polygons = []

    for pt in pts:
        k = (round(pt[0],4), round(pt[1],4))
        if k not in coord_to_region:
            polygons.append(Polygon())
            continue

        region_idx = coord_to_region[k][0]
        region_verts = vor.regions[region_idx]

        if -1 in region_verts: 
            polygons.append(Polygon())
            continue
        
        poly_coords = [vor.vertices[v] for v in region_verts]
        poly = Polygon(poly_coords)
        poly_clipped = poly.intersection(clip_box)
        polygons.append(poly_clipped)

    return polygons

def compute_area_per_lipid_voronoi(u, chol_sel, dopc_sel, folder, frame_stride=1):
    n_frames = len(u.trajectory[::frame_stride])
    results = np.zeros((n_frames, 4), dtype=np.float64)

    frame_idx = 0
    for ts in u.trajectory[::frame_stride]:
        lx, ly, lz, alpha, beta, gamma = ts.dimensions

        chol_coords = chol_sel.positions[:, :2]
        dopc_coords = dopc_sel.positions[:, :2]

        n_chol_res = len(chol_sel.residues)
        n_dopc_res = len(dopc_sel.residues)
        if n_chol_res + n_dopc_res == 0:
            results[frame_idx, :] = [ts.frame, 0.0, 0.0, 0.0]
            frame_idx += 1
            continue

        chol_type = np.zeros(chol_coords.shape[0], dtype=int)
        dopc_type = np.ones(dopc_coords.shape[0], dtype=int)

        all_coords = np.vstack((chol_coords, dopc_coords))
        all_types  = np.concatenate((chol_type, dopc_type))

        vor = voronoi_2d_pbc(all_coords, lx, ly, px=PX, py=PY)
        polygons = clip_voronoi_polygons_to_box(vor, all_coords, lx, ly, px=PX, py=PY)

        chol_areas = []
        dopc_areas = []
        for i, poly in enumerate(polygons):
            area_nm2 = poly.area if not poly.is_empty else 0.0
            if all_types[i] == 0:
                chol_areas.append(area_nm2)
            else:
                dopc_areas.append(area_nm2)

        if len(chol_areas) > 0:
            chol_area_per_lipid_nm2 = np.mean(chol_areas)
        else:
            chol_area_per_lipid_nm2 = 0.0

        if len(dopc_areas) > 0:
            dopc_area_per_lipid_nm2 = np.mean(dopc_areas)
        else:
            dopc_area_per_lipid_nm2 = 0.0

        chol_area_per_lipid_A2 = chol_area_per_lipid_nm2 * NM2_to_ANG2
        dopc_area_per_lipid_A2 = dopc_area_per_lipid_nm2 * NM2_to_ANG2

        n_chol = len(chol_areas)
        n_dopc = len(dopc_areas)
        total = n_chol + n_dopc
        if total > 0:
            avg_lipid_area = ((n_chol * chol_area_per_lipid_A2) + (n_dopc * dopc_area_per_lipid_A2)) / total
        else:
            avg_lipid_area = 0.0

        results[frame_idx, 0] = ts.frame
        results[frame_idx, 1] = chol_area_per_lipid_A2
        results[frame_idx, 2] = dopc_area_per_lipid_A2
        results[frame_idx, 3] = avg_lipid_area

        frame_idx += 1

    return results

def main():
    try:
        import shapely
    except ImportError:
        print("Shapely is required to clip Voronoi polygons! Install via 'pip install shapely'.")
        sys.exit(1)

    for folder in FOLDERS:
        path_gro = os.path.join(folder, GRO_FILE)
        path_xtc = os.path.join(folder, XTC_FILE)
        if not os.path.exists(path_gro) or not os.path.exists(path_xtc):
            print(f"Skipping folder {folder}: GRO or XTC not found.")
            continue

        print(f"\nProcessing folder: {folder}")

        u = mda.Universe(path_gro, path_xtc)

        chol_sel = u.select_atoms(f"resname {CHOL_RESNAME} and name {CHOL_HEAD_ATOM}")
        dopc_sel = u.select_atoms(f"resname {DOPC_RESNAME} and name {DOPC_HEAD_ATOM}")

        results = compute_area_per_lipid_voronoi(
            u, chol_sel, dopc_sel, folder, frame_stride=FRAME_STRIDE
        )

        outname = f"{folder}_area_voronoi.dat"
        header = "# Frame   CHOL_area(Å^2)   DOPC_area(Å^2)   Lipid_area(Å^2)"
        np.savetxt(outname, results, header=header)
        print(f"  -> Saved {outname}")

if __name__ == "__main__":
    main()

sns.set(style="whitegrid", font_scale=1, rc={"font.family": "Arial", "font.size": 18})

file_list = ["0_area_voronoi.dat", "10_area_voronoi.dat",
             "20_area_voronoi.dat", "30_area_voronoi.dat", "40_area_voronoi.dat"]
chol_percents = [0, 10, 20, 30, 40]

columns = [1, 2, 3]

avg_vals = {col: [] for col in columns}
std_vals = {col: [] for col in columns}

for filename in file_list:
    data = np.loadtxt(filename)
    for col in columns:
        avg_vals[col].append(np.mean(data[:, col]))
        std_vals[col].append(np.std(data[:, col]))

plot_titles = {
    1: "Average Cholesterol Area (Å²)",
    2: "Average DOPC Area (Å²)",
    3: "Average Lipid Area (Å²)"
}

for col in columns:
    plt.figure(figsize=(8, 6), dpi=300)
    plt.errorbar(chol_percents, avg_vals[col], yerr=std_vals[col],
                 fmt='o', capsize=5, markersize=8, color='b', ecolor='k', elinewidth=2, markeredgewidth=2)
    plt.xlabel("Cholesterol %")
    plt.ylabel("Average Area (Å²)")
    plt.title(plot_titles[col])
    plt.xticks(chol_percents)
    plt.tight_layout()
    plt.savefig(f'plot_{col}.png', dpi=300, bbox_inches='tight')
    plt.close()

