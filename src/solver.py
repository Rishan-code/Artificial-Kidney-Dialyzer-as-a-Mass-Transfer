# src/solver.py
import numpy as np
from scipy.optimize import fsolve
from parameters import *


# from mass_transfer import calculate_area, calculate_overall_Ko

def dialyzer_system(vars, Q_b, Q_d, C_b_in, C_d_in, Ko_A):
    """
    The system of non-linear equations to solve.
    vars is a list containing our unknowns: [N_total, C_b_out, C_d_out]
    """
    N, C_b_out, C_d_out = vars

    # Equation 1: Blood Mass Balance
    eq1 = N - Q_b * (C_b_in - C_b_out)

    # Equation 2: Dialysate Mass Balance
    eq2 = N - Q_d * (C_d_out - C_d_in)

    # Equation 3: LMCD Rate Equation
    delta_C1 = C_b_in - C_d_out
    delta_C2 = C_b_out - C_d_in

    # Prevent math domain errors if the differences become negative during solver iteration
    if delta_C1 <= 0 or delta_C2 <= 0 or delta_C1 == delta_C2:
        LMCD = 1e-6  # small placeholder
    else:
        LMCD = (delta_C1 - delta_C2) / np.log(delta_C1 / delta_C2)

    eq3 = N - Ko_A * LMCD

    return [eq1, eq2, eq3]


def analytical_concentration_profile(z, C_b_in, Ko_A, Q_b, Q_d):
    """
    Optional: Calculates the concentration at any point 'z' along the dialyzer.
    Solved using pure exponential terms (np.exp) to avoid hyperbolic instability.
    """
    # TODO: Implement the exponential decay formula for counter-current flow
    pass


def solve():
    """Main execution function."""
    # 1. Call your mass_transfer functions to get A and K_o
    # A = calculate_area(...)
    # Ko = calculate_overall_Ko(...)

    # For now, let's assume a placeholder value for Ko * A
    Ko_A_placeholder = 150 * (1e-6 / 60)

    # 2. Initial guesses for [N_total, C_b_out, C_d_out]
    initial_guesses = [1e-5, C_b_in * 0.5, C_b_in * 0.2]

    # 3. Run the solver
    solution = fsolve(dialyzer_system, initial_guesses,
                      args=(Q_b, Q_d, C_b_in, C_d_in, Ko_A_placeholder))

    print(f"Total Urea Removed (N): {solution[0]:.2e} mg/s")
    print(f"Blood Outlet Concentration: {solution[1]:.4f} mg/mL")
    print(f"Dialysate Outlet Concentration: {solution[2]:.4f} mg/mL")


if __name__ == "__main__":
    solve()