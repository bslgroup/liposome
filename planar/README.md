# Planar Bilayer Simulations

This directory contains simulation input files and analysis scripts for planar DOPC/cholesterol bilayers.

---

##  Description

Planar bilayers represent the **zero-curvature limit** and serve as a reference system for comparison with vesicle simulations.

---

## Contents



---

## Simulation Details

- MARTINI coarse-grained force field (v2.2)
- GROMACS (2024)
- Temperature: 310 K
- Time step: 20 fs
- Simulation time: ~7 µs

---

##  Analysis

Includes scripts for:
- Membrane thickness (global bilayer normal)
- Lipid interdigitation
- SASA
- Segmental order parameter (SCD)
- Lipid flip-flop

---

## Notes

- Planar systems exhibit uniform lipid packing
- Used as baseline for curvature comparison
- Large trajectory files are excluded
