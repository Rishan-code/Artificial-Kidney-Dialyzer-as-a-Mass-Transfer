# src/mass_transfer.py
import numpy as np
from parameters import *


def calculate_area():
    """Calculates the total membrane surface area."""
    return n_fibers * np.pi * d_i * L


def calculate_diffusivities():
    """
    Uses the Stokes-Einstein equation to calculate free-water diffusivity,
    and structural parameters to calculate effective membrane diffusivity.
    """
    # 1. Stokes-Einstein for Diffusivity in Water/Plasma (m^2/s)
    D_water = (kb_boltzmann * T_kelvin) / (6 * np.pi * mu_water * r_s)

    # 2. Membrane Diffusivity based on polymer structure
    D_m = D_water * (porosity / tortuosity)

    return D_water, D_m


def calculate_kb(Q_b_m3_s, D_blood):
    """
    Calculates the blood-side mass transfer coefficient (k_b) using the
    Leveque correlation for laminar flow inside hollow fibers.
    """
    # Dynamic viscosity of blood based on Hematocrit
    mu_blood = mu_water * (1 + 2.5 * Hct)
    rho_blood = 1060.0  # kg/m^3

    # Cross-sectional area of all fibers
    A_cross = n_fibers * np.pi * (d_i / 2) ** 2
    velocity = Q_b_m3_s / A_cross

    # Dimensionless numbers
    Re = (rho_blood * velocity * d_i) / mu_blood
    Sc = mu_blood / (rho_blood * D_blood)

    # Leveque solution for Sherwood number in developing laminar pipe flow
    # Graetz number (Gz) = Re * Sc * (d_i / L)
    Gz = Re * Sc * (d_i / L)

    if Gz > 10:
        Sh = 1.62 * (Gz ** (1 / 3))
    else:
        Sh = 3.66  # Asymptotic limit for fully developed flow

    k_b = (Sh * D_blood) / d_i
    return k_b


def calculate_kd(Q_d_m3_s, D_dialysate):
    """
    Calculates the dialysate-side mass transfer coefficient (k_d).
    Uses a simplified empirical correlation for shell-side flow.
    """
    # Assume dialysate properties are close to water
    rho_d = 993.0
    mu_d = mu_water

    # Superficial velocity on the shell side (simplified assumption)
    # Assuming a shell diameter roughly wrapping the fiber bundle
    D_shell = np.sqrt(n_fibers) * d_i * 1.5
    A_shell = (np.pi / 4) * (D_shell ** 2) - (n_fibers * np.pi * (d_i / 2) ** 2)
    velocity_d = Q_d_m3_s / A_shell

    Re_d = (rho_d * velocity_d * d_i) / mu_d
    Sc_d = mu_d / (rho_d * D_dialysate)

    # Typical empirical Sherwood correlation for cross/shell flow
    Sh_d = 0.5 * (Re_d ** 0.6) * (Sc_d ** 0.33)

    k_d = (Sh_d * D_dialysate) / d_i
    return k_d


def calculate_overall_Ko(Q_b_m3_s, Q_d_m3_s):
    """
    Master function to calculate the overall mass transfer coefficient (K_o).
    """
    D_water, D_m = calculate_diffusivities()

    k_b = calculate_kb(Q_b_m3_s, D_water)
    k_d = calculate_kd(Q_d_m3_s, D_water)  # Assuming dialysate D is similar to water

    # Resistance-in-series
    R_blood = 1 / k_b
    R_membrane = delta / D_m
    R_dialysate = 1 / k_d

    R_total = R_blood + R_membrane + R_dialysate
    K_o = 1 / R_total

    return K_o, R_blood, R_membrane, R_dialysate