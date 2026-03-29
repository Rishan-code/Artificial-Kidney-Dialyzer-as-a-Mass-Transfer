# src/solver.py
import numpy as np
from scipy.optimize import fsolve
from parameters import *
from mass_transfer import calculate_area, calculate_overall_Ko  # <-- NEW IMPORT


def dialyzer_system(vars, Q_b, Q_d, C_b_in, C_d_in, Ko_A):
    # ... (Keep this function exactly as you have it) ...
    N, C_b_out, C_d_out = vars
    eq1 = N - Q_b * (C_b_in - C_b_out)
    eq2 = N - Q_d * (C_d_out - C_d_in)

    delta_C1 = C_b_in - C_d_out
    delta_C2 = C_b_out - C_d_in
    if delta_C1 <= 0 or delta_C2 <= 0 or delta_C1 == delta_C2:
        LMCD = 1e-6
    else:
        LMCD = (delta_C1 - delta_C2) / np.log(delta_C1 / delta_C2)
    eq3 = N - Ko_A * LMCD
    return [eq1, eq2, eq3]


def solve():
    """Main execution function."""
    # 1. Calculate Area (Constant for the device)
    A = calculate_area()

    # 2. Calculate dynamic Ko based on baseline flow rates in parameters.py
    Ko, R_blood, R_membrane, R_dialysate = calculate_overall_Ko(Q_b, Q_d)

    # Multiply them to get the real overall capacity
    Ko_A_actual = Ko * A

    print(f"--- Resistance Breakdown ---")
    print(
        f"Blood Film: {R_blood / (1 / Ko) * 100:.1f}% | Membrane: {R_membrane / (1 / Ko) * 100:.1f}% | Dialysate: {R_dialysate / (1 / Ko) * 100:.1f}%")
    print("-" * 30)

    # 3. Initial guesses and Solver
    initial_guesses = [1e-5, C_b_in * 0.5, C_b_in * 0.2]
    solution = fsolve(dialyzer_system, initial_guesses,
                      args=(Q_b, Q_d, C_b_in, C_d_in, Ko_A_actual))

    print(f"Total Urea Removed (N): {solution[0]:.2e} mg/s")
    print(f"Blood Outlet Concentration: {solution[1]:.4f} mg/mL")
    print(f"Dialysate Outlet Concentration: {solution[2]:.4f} mg/mL")


if __name__ == "__main__":
    solve()