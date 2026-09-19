# Adaptive Optics — Deformable Mirror Design Optimization

Machine learning-driven parameter optimization for a Shape Memory Alloy (SMA)-actuated deformable mirror, aimed at improving wavefront correction accuracy for adaptive optics applications.

![Status](https://img.shields.io/badge/status-complete-brightgreen) ![MATLAB](https://img.shields.io/badge/MATLAB-simulation-orange) ![Python](https://img.shields.io/badge/python-scikit--learn-blue) ![License](https://img.shields.io/badge/license-MIT-lightgrey)

---

## Overview

A deformable mirror corrects optical wavefront aberrations (defocus, astigmatism, higher-order Zernike modes) by physically reshaping its surface via actuators. This project models the relationship between a mirror's geometric design parameters and its wavefront-correction accuracy, then uses machine learning to search that parameter space for a design that outperforms the original hand-picked configuration — and to solve the inverse problem: given a target correction accuracy, what geometry produces it.

## Design Parameters

The mirror's mechanical geometry is defined by 5 parameters, each varied across a discrete range for simulation:

| Parameter | Description | Range Tested |
|---|---|---|
| Arm Length | Length of each actuator arm | 37.5 – 42.5 |
| Arm Angle | Angular spacing between arms | 7 – 10.5 |
| Inner Diameter | Diameter of the central deforming region | 55 – 60 |
| Thickness | Mirror substrate thickness | 0.201 mm (fixed) |
| OACD | Actuator Contact Diameter | 4 – 6 |
| Fillet Radius | Corner radius at arm junctions | 1 – 5 |

## Pipeline

```
MATLAB: Parametric mirror geometry + influence matrix generation
   → Zernike-mode wavefront simulation (defocus, astig45, astig0, higher-order Z)
   → Python: Correlation analysis of design parameters vs. correction accuracy
   → Regression modeling (Elastic-Net, Ridge, Polynomial, Random Forest)
   → Inverse optimization (scipy.minimize) — reverse-engineer optimal geometry
   → Validated final design vs. baseline
```

## What Was Done

**1. Simulation & Influence Matrix (MATLAB)**
Built a parametric model of the mirror geometry and generated an influence matrix mapping actuator inputs to surface deformation, used to simulate wavefront correction across multiple Zernike modes.

**2. Correlation Analysis**
Analyzed relationships between design parameters and correction accuracy:
- Inner diameter and arm length showed strong negative correlation with several output modes
- OACD and arm angle showed moderate positive correlation
- Fillet radius showed weak correlation, indicating low design sensitivity
- Several output modes (e.g. astig45, astig0, and related Zernike terms) were highly correlated with each other, suggesting they could be modeled jointly

**3. Regression Modeling**
Trained and compared multiple models to predict correction accuracy from geometry parameters, evaluated on bias/variance trade-off:

| Model | Bias | Variance |
|---|---|---|
| Random Forest | 8.93 | 59.60 |
| Polynomial RR (deg 3) | 17.17 | 106.49 |
| Polynomial RR (deg 2) | 26.66 | 66.13 |
| Elastic-Net | 41.38 | 72.63 |
| Linear Regression | 57.33 | 89.23 |
| Ridge Regression | 58.27 | 91.59 |

Random Forest gave the best bias-variance balance and was used for final parameter selection.

**4. Inverse Optimization**
Used `scipy.minimize` to solve the inverse design problem — given a target correction accuracy (Y_target), find the geometry parameters (X) that best achieve it, by minimizing squared error between model prediction and target across a bounded parameter search.

## Results

Comparing the original hand-tuned design against the ML-optimized geometry:

| Zernike Mode | Initial Model | ML-Optimized |
|---|---|---|
| Defocus | 68.0 | **82.58** |
| Astig 45° | 72.9 | **86.30** |
| Astig 0° | 71.0 | **86.98** |
| Z = 17 | 65.6 | **90.31** |

**Final optimized geometry:** Arm Length 37.5, Arm Angle 10, Inner Diameter 30, OACD 6, Fillet Radius 1 — a measurable improvement in wavefront-correction accuracy across all tested modes over the original design.

## Repository Structure

```
├── matlab/                # Mirror geometry, influence matrix, wavefront simulation
├── python/
│   ├── correlation_analysis.py
│   ├── regression_models.py     # Elastic-Net, Ridge, Polynomial, Random Forest
│   └── inverse_optimization.py  # scipy.minimize-based inverse design
├── results/                # Correlation matrices, model comparison, wavefront plots
└── README.md
```

## Tech Stack

MATLAB (mirror simulation, wavefront modeling) · Python (Pandas, Scikit-learn, SciPy) · PTC Creo (mechanical design)

## Future Work

- Extend the parameter search to include material properties (SMA composition, Kapton layer thickness)
- Validate ML-optimized geometry against physical prototype measurements
- Explore multi-output regression to jointly predict correlated Zernike modes

## License

MIT
