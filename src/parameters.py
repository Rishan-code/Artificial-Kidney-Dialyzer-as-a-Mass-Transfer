# src/parameters.py

# --- Geometric Parameters of the Dialyzer ---
d_i = 200e-6        # Inner diameter of a single hollow fiber (meters)
delta = 40e-6       # Membrane thickness (meters)
L = 0.25            # Active length of the dialyzer (meters)
n_fibers = 10000    # Total number of hollow fibers

# --- Operating Conditions ---
Q_b_ml_min = 300    # Blood flow rate (mL/min)
Q_d_ml_min = 500    # Dialysate flow rate (mL/min)

# Convert flow rates to standard SI units (m^3/s) for calculations
Q_b = Q_b_ml_min * (1e-6 / 60)
Q_d = Q_d_ml_min * (1e-6 / 60)

# --- Inlet Concentrations ---
C_b_in = 1.5        # Inlet urea concentration in blood (mg/mL)
C_d_in = 0.0        # Clean dialysate entering device (mg/mL)

# --- Diffusivities (m^2/s) - Placeholders to be researched ---
D_blood = 1.0e-9    # Diffusivity of urea in blood
D_m = 0.5e-9        # Effective diffusivity of urea in the membrane
D_dialysate = 1.2e-9 # Diffusivity of urea in dialysate