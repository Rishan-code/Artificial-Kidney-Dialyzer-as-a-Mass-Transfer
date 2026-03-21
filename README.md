# Hemodialyzer Mass Transfer Simulation

## Overview
This repository contains a computational model of an artificial kidney (hemodialyzer) developed for a Separation Processes course project. The model simulates the removal of metabolic waste (urea) from blood across a semi-permeable hollow-fiber membrane using steady-state mass transfer principles.

The simulation acts as a macroscopic process model, treating the dialyzer analogously to a counter-current shell-and-tube heat exchanger, utilizing the Log Mean Concentration Difference (LMCD) and the resistance-in-series model to predict urea clearance.

## Core Methodology
The fundamental physics governing this simulation include:
* **Fick's Law of Diffusion:** For basic transport across the membrane.
* **Resistance-in-Series Model:** Calculating the overall mass transfer coefficient ($K_o$) by combining blood-side film resistance ($k_b$), membrane resistance, and dialysate-side film resistance ($k_d$).
* **LMCD Approach:** Solving the non-linear concentration gradients inherent in counter-current flow operations.

## Repository Structure

├── src/
│   ├── __init__.py
│   ├── parameters.py      # Standard physiological and geometric constants
│   ├── mass_transfer.py   # Fluid dynamics and Sherwood number correlations
│   └── solver.py          # Numerical solver for simultaneous mass balances
├── notebooks/
│   └── sensitivity_analysis.ipynb  # Visualization of flow rates vs. clearance
├── requirements.txt
└── README.md


## Getting Started

### Prerequisites
Ensure you have Python 3.8+ installed. The primary dependencies are `numpy`, `scipy`, and `matplotlib` (for the Jupyter notebooks).

### Installation
1. Clone the repository:
   ```bash
   git clone [https://github.com/Rishan-code/dialyzer-simulation.git](https://github.com/Rishan-code/dialyzer-simulation.git)