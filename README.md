# Liposome Curvature–Cholesterol MD Reproducibility

This repository contains input files and analysis scripts for reproducing coarse-grained molecular dynamics (CG-MD) simulations of DOPC/cholesterol membranes.

---

##  Scientific Objective

This study investigates how **cholesterol concentration** and **membrane curvature** affect the structural and dynamic properties of lipid bilayers.

Two systems are considered:
- Planar bilayers (flat membranes)
- Spherical vesicles (curved membranes)

---

## Repository Structure



---

## Simulation Details

- **Force field:** MARTINI coarse-grained (v2.2)
- **Software:** GROMACS (2024)
- **Temperature:** 310 K
- **Time step:** 20 fs
- **Simulation time:** ~7 microseconds

### Workflow:
1. Energy minimization
2. NVT equilibration
3. NPT equilibration
4. Production MD

---

##  Analysis

Analysis performed using:
- Python (MDAnalysis, NumPy, Matplotlib)
- GROMACS tools (`gmx sasa`)

Computed properties:
- Membrane thickness
- Lipid interdigitation
- Solvent-accessible surface area (SASA)
- Segmental order parameter (SCD)
- Lipid flip-flop dynamics
- Cholesterol distribution

---

## Notes

- Large trajectory files (`.xtc`, `.trr`) are not included
- Only input files and scripts required for reproducibility are provided
- MARTINI simulations represent accelerated dynamics (relative trends)

---

##  Citation

Khodadadi et al. 
*Interplay Between Cholesterol Concentration and Membrane Curvature in Liposomes Revealed by Molecular Dynamics Simulations*

---

## Authors

Ehsan Khodadadi 
Mortaza Derakhshani-Molayousefi 
Ehsaneh Khodadadi 
Mahmoud Moradi 

University of Arkansas

---

## Contact

moradi@uark.edu
