# src/parameters.py
import numpy as np

# --- 1. Solute Properties ---
# Dictionary of solutes to simulate: Name -> {'Mw': g/mol, 'r_s': Stokes radius in m}
SOLUTES = {
    'Urea': {'Mw': 60.06, 'r_s': 2.6e-10},
    'VitB12': {'Mw': 1355.0, 'r_s': 8.5e-10}, # Middle molecule
    'Albumin': {'Mw': 66000.0, 'r_s': 35.8e-10} # Large protein
}

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
Q_uf_ml_min = 10.0      # Ultrafiltration rate (water removal)

Q_b = Q_b_ml_min * (1e-6 / 60)
Q_d = Q_d_ml_min * (1e-6 / 60)
Q_uf = Q_uf_ml_min * (1e-6 / 60)

# Inlet concentrations (mg/mL)
C_b_in = {
    'Urea': 1.5,
    'VitB12': 0.02,
    'Albumin': 40.0
}
C_d_in = 0.0

# --- Base physical constants for calculations ---
kb_boltzmann = 1.38e-23 # Boltzmann constant
mu_water = 6.9e-4       # Dynamic viscosity of water at 37C (Pa*s)

# Example: Dynamic calculation of Blood Viscosity based on Hematocrit
mu_blood = mu_water * (1 + 2.5 * Hct)