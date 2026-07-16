# Vesicle (Spherical Bilayer) Simulations

This directory contains simulation input files and analysis scripts for spherical liposome systems.

---

##  Description

Vesicle systems introduce **membrane curvature**, which affects:
- lipid packing
- interleaflet interactions
- membrane dynamics

These systems are used to study curvature-dependent membrane behavior.

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
- Membrane thickness (radial definition)
- Lipid interdigitation (curvature-dependent)
- SASA
- Segmental order parameter (SCD)
- Lipid flip-flop dynamics
- Vesicle structural properties

---

## Notes

- Curvature introduces structural heterogeneity
- Results represent global averages over vesicle surface
- Used to compare with planar bilayers
- Large trajectory files are excluded
