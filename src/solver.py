# src/solver.py
import numpy as np
from scipy.integrate import solve_bvp
import sys
import os

# Add the src structure so it can be run standalone
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from parameters import *
try:
    from mass_transfer import calculate_area, calculate_overall_Ko, calculate_sieving_coefficient
except ImportError:
    from .mass_transfer import calculate_area, calculate_overall_Ko, calculate_sieving_coefficient

def solve_solute_profile(solute_name, qb_in, qd_in, quf):
    """
    Solves the 1D spatial concentration profile for a specific solute
    using a Boundary Value Problem (BVP) solver.
    """
    solute_params = SOLUTES[solute_name]
    r_s = solute_params['r_s']
    cb_in_val = C_b_in[solute_name]
    
    A = calculate_area()
    S = calculate_sieving_coefficient(r_s)
    
    # We will use average flows to calculate a constant Ko, or recalculate Ko along length.
    # For a macroscopic model, constant Ko based on mean flow is standard.
    Ko, _, _, _ = calculate_overall_Ko(qb_in - quf/2, qd_in + quf/2, r_s)
    
    def ode_system(z, y):
        # y[0] = C_b, y[1] = C_d
        Cb = y[0]
        Cd = y[1]
        
        # Flows at position z
        # Blood flows z=0 to z=L
        Qb_z = qb_in - quf * (z / L)
        # Dialysate flows z=L to z=0, so its positive magnitude at z is:
        Qd_z = qd_in + quf * ((L - z) / L)
        
        dCb_dz = (1 / Qb_z) * (-Ko * (A / L) * (Cb - Cd) + (quf / L) * (1 - S) * Cb)
        dCd_dz = (1 / Qd_z) * (-Ko * (A / L) * (Cb - Cd) - (quf / L) * (S * Cb - Cd))
        
        return np.vstack((dCb_dz, dCd_dz))
        
    def boundary_conditions(ya, yb):
        # ya corresponds to z=0. Blood enters here.
        # yb corresponds to z=L. Dialysate enters here.
        return np.array([
            ya[0] - cb_in_val,  # Cb(0) = Cb_in
            yb[1] - C_d_in      # Cd(L) = Cd_in (= 0)
        ])
        
    # Initial mesh
    z_mesh = np.linspace(0, L, 50)
    # Initial guess: linear drop for blood, linear rise for dialysate
    y_guess = np.zeros((2, z_mesh.size))
    y_guess[0] = np.linspace(cb_in_val, cb_in_val*0.2, z_mesh.size)
    y_guess[1] = np.linspace(cb_in_val*0.5, 0, z_mesh.size)
    
    sol = solve_bvp(ode_system, boundary_conditions, z_mesh, y_guess)
    
    if not sol.success:
        print(f"Warning: Convergence failed for {solute_name}")
        
    # Calculate Clearance (K)
    # K = [Qb_in * Cb_in - Qb_out * Cb_out] / Cb_in
    Cb_out = sol.y[0, -1]
    Qb_out = qb_in - quf
    mass_in = qb_in * cb_in_val
    mass_out = Qb_out * Cb_out
    mass_removed = mass_in - mass_out
    
    K_m3_s = mass_removed / cb_in_val if cb_in_val > 0 else 0
    K_ml_min = K_m3_s * 60 * 1e6
    
    return {
        'z': sol.x,
        'Cb_profile': sol.y[0],
        'Cd_profile': sol.y[1],
        'Clearance_ml_min': K_ml_min,
        'Sieving_Coefficient': S,
        'Cb_out': Cb_out
    }

def solve_all(print_results=True):
    """Solves the profiles for all solutes and prints results."""
    results = {}
    
    if print_results:
        print("="*50)
        print("HEMODIALYZER SIMULATION RESULTS")
        print("="*50)
        print(f"Operating Conditions: Qb = {Q_b_ml_min} mL/min | Qd = {Q_d_ml_min} mL/min | Quf = {Q_uf_ml_min} mL/min")
        
    for solute in SOLUTES.keys():
        res = solve_solute_profile(solute, Q_b, Q_d, Q_uf)
        results[solute] = res
        
        if print_results:
            print(f"\n--- {solute} ---")
            print(f"Sieving Coefficient (S): {res['Sieving_Coefficient']:.4f}")
            print(f"Clearance (K):           {res['Clearance_ml_min']:.1f} mL/min")
            print(f"Outlet Blood Conc.:      {res['Cb_out']:.4f} mg/mL (from {C_b_in[solute]} mg/mL)")
            
    return results

if __name__ == "__main__":
    solve_all()