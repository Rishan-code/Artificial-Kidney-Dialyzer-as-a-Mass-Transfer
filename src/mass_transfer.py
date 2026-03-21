# src/mass_transfer.py
import numpy as np
from parameters import *

def calculate_area(d_i, L, n_fibers):
    """Calculates the total membrane surface area."""
    # TODO: Implement formula (n * pi * d_i * L)
    pass

def calculate_kb(Q_b, d_i, n_fibers, D_blood):
    """
    Calculates the blood-side mass transfer coefficient (k_b).
    Requires calculating blood velocity, Reynolds number, and Sherwood number.
    """
    # TODO: Add fluid dynamics correlations here
    pass

def calculate_kd(Q_d, D_dialysate):
    """
    Calculates the dialysate-side mass transfer coefficient (k_d) on the shell side.
    """
    # TODO: Add shell-side Sherwood number correlation here
    pass

def calculate_overall_Ko(k_b, k_d, D_m, delta):
    """
    Calculates the overall mass transfer coefficient (K_o) using the
    resistance-in-series model.
    """
    # TODO: Implement 1/Ko = 1/kb + delta/Dm + 1/kd
    pass