# src/parameters.py
import numpy as np

# --- 1. Solute Properties (Variable Solute: Urea) ---
Mw = 60.06              # Molecular weight (g/mol)
r_s = 2.6e-10           # Stokes radius of the molecule (meters)
protein_bound = 0.0     # Fraction of solute bound to large proteins (0 to 1)

# --- 2. Polymer Membrane Structural Variables ---
d_i = 200e-6            # Inner diameter of fiber (m)
delta = 40e-6           # Membrane thickness (m)
L = 0.25                # Length (m)
n_fibers = 10000        # Number of fibers
porosity = 0.75         # Membrane porosity (epsilon)
tortuosity = 2.5        # Membrane tortuosity (tau)
pore_radius = 3.0e-9    # Average pore radius (meters)

# --- 3. Fluid & Patient Variables ---
T_celsius = 37.0        # Body temperature (C)
T_kelvin = T_celsius + 273.15
Hct = 0.40              # Hematocrit (40% red blood cells)

# Flow rates (converted to m^3/s)
Q_b_ml_min = 300
Q_d_ml_min = 500
Q_b = Q_b_ml_min * (1e-6 / 60)
Q_d = Q_d_ml_min * (1e-6 / 60)

# Inlet concentrations (mg/mL)
C_b_in = 1.5
C_d_in = 0.0

# --- Base physical constants for calculations ---
kb_boltzmann = 1.38e-23 # Boltzmann constant
mu_water = 6.9e-4       # Dynamic viscosity of water at 37C (Pa*s)

# Example: Dynamic calculation of Blood Viscosity based on Hematocrit
mu_blood = mu_water * (1 + 2.5 * Hct)