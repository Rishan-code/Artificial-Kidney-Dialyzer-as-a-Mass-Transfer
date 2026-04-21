#!/usr/bin/env python3
"""
demo.py — Hemodialyzer Process Simulator Demo
=============================================
Run this script for the Apr 22 class demonstration.

Usage:
    python demo.py                          # default parameters
    python demo.py --qb 250 --qd 600       # custom blood/dialysate flow
    python demo.py --interactive            # interactive mode (tweak params live)
"""

import argparse
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

# ── Setup paths ──────────────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from parameters import SOLUTES, C_b_in, C_d_in
from solver import solve_solute_profile
from mass_transfer import (
    calculate_overall_Ko,
    calculate_area,
    calculate_diffusivities,
    calculate_sieving_coefficient,
    calculate_kb,
    calculate_kd,
)

# ── Color palette ────────────────────────────────────────────────────────
COLORS = {
    "Urea":    {"blood": "#E63946", "dial": "#457B9D", "op": "#E63946"},
    "VitB12":  {"blood": "#2A9D8F", "dial": "#264653", "op": "#2A9D8F"},
    "Albumin": {"blood": "#E9C46A", "dial": "#F4A261", "op": "#E9C46A"},
}
BG_COLOR     = "#0F1117"
PANEL_COLOR  = "#1A1B26"
GRID_COLOR   = "#2A2B3D"
TEXT_COLOR   = "#E8E8F0"
ACCENT       = "#7AA2F7"
EQUIL_COLOR  = "#FF6B6B"


# ══════════════════════════════════════════════════════════════════════════
#  CORE SOLVER WRAPPER
# ══════════════════════════════════════════════════════════════════════════
def run_simulation(qb_ml_min, qd_ml_min, quf_ml_min=10.0):
    """Run the BVP solver for all solutes at the given operating conditions."""
    qb  = qb_ml_min  * (1e-6 / 60)
    qd  = qd_ml_min  * (1e-6 / 60)
    quf = quf_ml_min * (1e-6 / 60)

    results = {}
    for solute in SOLUTES:
        res = solve_solute_profile(solute, qb, qd, quf)
        results[solute] = res
    return results


# ══════════════════════════════════════════════════════════════════════════
#  RESISTANCE BREAKDOWN HELPER
# ══════════════════════════════════════════════════════════════════════════
def get_resistance_breakdown(qb_ml_min, qd_ml_min):
    """Return individual resistances for each solute."""
    qb = qb_ml_min * (1e-6 / 60)
    qd = qd_ml_min * (1e-6 / 60)
    breakdown = {}
    for name, props in SOLUTES.items():
        Ko, Rb, Rm, Rd = calculate_overall_Ko(qb, qd, props['r_s'])
        Rtot = Rb + Rm + Rd
        breakdown[name] = {
            "Ko": Ko,
            "R_blood_%":     100 * Rb / Rtot,
            "R_membrane_%":  100 * Rm / Rtot,
            "R_dialysate_%": 100 * Rd / Rtot,
        }
    return breakdown


# ══════════════════════════════════════════════════════════════════════════
#  PRETTY TERMINAL OUTPUT
# ══════════════════════════════════════════════════════════════════════════
def print_results(results, qb, qd, quf):
    """Print a clean summary table to the terminal."""
    hline = "=" * 70
    print(f"\n{hline}")
    print("  HEMODIALYZER PROCESS SIMULATOR - RESULTS")
    print(hline)
    print(f"  Operating Conditions:")
    print(f"    Blood flow  (Qb)  = {qb:>6.0f} mL/min")
    print(f"    Dialysate   (Qd)  = {qd:>6.0f} mL/min")
    print(f"    Ultrafilt.  (Quf) = {quf:>6.1f} mL/min")
    print(hline)
    print(f"  {'Solute':<12} {'MW (Da)':>10} {'Sieving':>10} {'Clearance':>14} {'Cb_out':>12}")
    print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*14} {'-'*12}")
    
    for name, res in results.items():
        mw = SOLUTES[name]['Mw']
        print(
            f"  {name:<12} {mw:>10.0f} {res['Sieving_Coefficient']:>10.4f}"
            f" {res['Clearance_ml_min']:>11.1f} mL/min"
            f" {res['Cb_out']:>9.4f} mg/mL"
        )
    print(hline)

    # Resistance breakdown
    bd = get_resistance_breakdown(qb, qd)
    print(f"\n  RESISTANCE BREAKDOWN (% of total)")
    print(f"  {'Solute':<12} {'Ko (m/s)':>12} {'Blood %':>10} {'Membrane %':>12} {'Dialysate %':>13}")
    print(f"  {'-'*12} {'-'*12} {'-'*10} {'-'*12} {'-'*13}")
    for name, info in bd.items():
        print(
            f"  {name:<12} {info['Ko']:>12.2e}"
            f" {info['R_blood_%']:>9.1f}%"
            f" {info['R_membrane_%']:>11.1f}%"
            f" {info['R_dialysate_%']:>12.1f}%"
        )
    print(hline + "\n")


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 1 — CONCENTRATION PROFILES ALONG FIBER LENGTH
# ══════════════════════════════════════════════════════════════════════════
def plot_concentration_profiles(results, ax_list):
    """Plot Cb(z) and Cd(z) for each solute on separate subplots."""
    for ax, (name, res) in zip(ax_list, results.items()):
        z_cm = res["z"] * 100  # convert m → cm
        cb = res["Cb_profile"]
        cd = res["Cd_profile"]

        ax.fill_between(z_cm, cb, alpha=0.15, color=COLORS[name]["blood"])
        ax.fill_between(z_cm, cd, alpha=0.10, color=COLORS[name]["dial"])
        ax.plot(z_cm, cb, lw=2.5, color=COLORS[name]["blood"], label="$C_b(z)$ — Blood")
        ax.plot(z_cm, cd, lw=2.5, color=COLORS[name]["dial"],  label="$C_d(z)$ — Dialysate",
                linestyle="--")

        ax.set_title(name, fontsize=13, fontweight="bold", color=TEXT_COLOR, pad=8)
        ax.set_xlabel("Fiber length z (cm)", fontsize=10, color=TEXT_COLOR)
        ax.set_ylabel("Concentration (mg/mL)", fontsize=10, color=TEXT_COLOR)
        ax.legend(fontsize=8, loc="best", framealpha=0.6,
                  facecolor=PANEL_COLOR, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR)
        ax.set_facecolor(PANEL_COLOR)
        ax.tick_params(colors=TEXT_COLOR, labelsize=9)
        ax.grid(True, color=GRID_COLOR, linewidth=0.5, alpha=0.5)
        for spine in ax.spines.values():
            spine.set_color(GRID_COLOR)

        # Annotate clearance
        K = res["Clearance_ml_min"]
        ax.text(0.97, 0.95, f"K = {K:.1f} mL/min",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=10, fontweight="bold", color=ACCENT,
                bbox=dict(boxstyle="round,pad=0.3", facecolor=BG_COLOR,
                          edgecolor=ACCENT, alpha=0.85))


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 2 — OPERATING DIAGRAM  (Cb vs Cd)  +  EQUILIBRIUM LINE
# ══════════════════════════════════════════════════════════════════════════
def plot_operating_diagram(results, ax):
    """
    The dialysis analogy of a McCabe–Thiele diagram.
    
    • Equilibrium line: Cb = Cd  (45° line — no driving force)
    • Operating curves: parametric (Cd(z), Cb(z)) for each solute
    • The horizontal gap between operating curve and equilibrium line 
      represents the local mass-transfer driving force (Cb − Cd).
    """
    max_conc = 0
    handles = []

    for name, res in results.items():
        cb = res["Cb_profile"]
        cd = res["Cd_profile"]
        color = COLORS[name]["op"]

        # Skip Albumin if its clearance is essentially zero (clutters the plot)
        if name == "Albumin" and res["Clearance_ml_min"] < 0.5:
            continue

        max_conc = max(max_conc, cb.max(), cd.max())

        # Operating curve
        ax.plot(cd, cb, lw=2.5, color=color, zorder=3)
        
        # Mark inlet and outlet
        # Blood inlet: z=0 → first point (high Cb, low Cd)
        ax.scatter(cd[0],  cb[0],  s=80, color=color, zorder=5, edgecolors="white", linewidths=1.2)
        # Blood outlet: z=L → last point (low Cb, high Cd)
        ax.scatter(cd[-1], cb[-1], s=80, color=color, zorder=5, edgecolors="white", linewidths=1.2,
                   marker="s")

        # Annotate the driving force at the midpoint
        mid = len(cb) // 2
        cb_mid, cd_mid = cb[mid], cd[mid]
        driving_force = cb_mid - cd_mid
        ax.annotate(
            "", xy=(cb_mid, cb_mid), xytext=(cd_mid, cb_mid),
            arrowprops=dict(arrowstyle="<->", color=color, lw=1.5, alpha=0.7),
        )
        ax.text(
            (cd_mid + cb_mid) / 2, cb_mid + max_conc * 0.02,
            f"dC = {driving_force:.3f}",
            ha="center", fontsize=8, color=color, fontweight="bold",
        )

        handles.append(Line2D([0], [0], color=color, lw=2.5, label=name))

    # ── Equilibrium line (45°) ──
    lim = max_conc * 1.1
    ax.plot([0, lim], [0, lim], ls="--", lw=2, color=EQUIL_COLOR, alpha=0.8, zorder=2,
            label="Equilibrium ($C_b = C_d$)")
    handles.append(Line2D([0], [0], color=EQUIL_COLOR, lw=2, ls="--",
                          label="Equilibrium ($C_b = C_d$)"))

    # ── Shade the region between operating curve and equilibrium for Urea ──
    if "Urea" in results:
        cb_u = results["Urea"]["Cb_profile"]
        cd_u = results["Urea"]["Cd_profile"]
        # Sort by Cd for clean fill
        sort_idx = np.argsort(cd_u)
        cd_sorted = cd_u[sort_idx]
        cb_sorted = cb_u[sort_idx]
        ax.fill_betweenx(cb_sorted, cd_sorted, cb_sorted,
                         alpha=0.08, color=COLORS["Urea"]["op"],
                         label="_driving force region")

    ax.set_xlim(left=-0.01)
    ax.set_ylim(bottom=-0.01)
    ax.set_xlabel("$C_d$ — Dialysate concentration (mg/mL)", fontsize=11, color=TEXT_COLOR)
    ax.set_ylabel("$C_b$ — Blood concentration (mg/mL)", fontsize=11, color=TEXT_COLOR)
    ax.set_title("Operating Diagram — Cb vs Cd\n(Equilibrium = 45° line, gap = driving force)",
                 fontsize=13, fontweight="bold", color=TEXT_COLOR, pad=10)
    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.7,
              facecolor=PANEL_COLOR, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR)
    ax.set_facecolor(PANEL_COLOR)
    ax.tick_params(colors=TEXT_COLOR, labelsize=9)
    ax.grid(True, color=GRID_COLOR, linewidth=0.5, alpha=0.5)
    ax.set_aspect("equal", adjustable="datalim")
    for spine in ax.spines.values():
        spine.set_color(GRID_COLOR)


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 3 — RESISTANCE PIE CHARTS
# ══════════════════════════════════════════════════════════════════════════
def plot_resistance_pies(qb, qd, ax_list):
    """Pie chart showing resistance breakdown for each solute."""
    bd = get_resistance_breakdown(qb, qd)
    pie_colors = ["#E63946", "#457B9D", "#2A9D8F"]
    labels = ["Blood film", "Membrane", "Dialysate film"]

    for ax, (name, info) in zip(ax_list, bd.items()):
        sizes = [info["R_blood_%"], info["R_membrane_%"], info["R_dialysate_%"]]
        wedges, texts, autotexts = ax.pie(
            sizes, labels=labels, autopct="%1.1f%%", startangle=140,
            colors=pie_colors, textprops={"fontsize": 8, "color": TEXT_COLOR},
            wedgeprops={"edgecolor": BG_COLOR, "linewidth": 1.5},
            pctdistance=0.6,
        )
        for t in autotexts:
            t.set_fontsize(8)
            t.set_fontweight("bold")
        ax.set_title(f"{name}\n$K_o$ = {info['Ko']:.2e} m/s",
                     fontsize=11, fontweight="bold", color=TEXT_COLOR, pad=6)
        ax.set_facecolor(PANEL_COLOR)


# ══════════════════════════════════════════════════════════════════════════
#  MASTER FIGURE BUILDER
# ══════════════════════════════════════════════════════════════════════════
def build_demo_figure(results, qb, qd, quf):
    """Assemble the full 3-row demo figure."""
    fig = plt.figure(figsize=(18, 18), facecolor=BG_COLOR)
    fig.suptitle(
        f"Hemodialyzer Process Simulator   |   $Q_b$={qb:.0f}   $Q_d$={qd:.0f}   $Q_{{uf}}$={quf:.0f}  mL/min",
        fontsize=16, fontweight="bold", color=TEXT_COLOR, y=0.98,
    )

    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35,
                  left=0.06, right=0.96, top=0.93, bottom=0.04)

    # Row 1 — Concentration profiles (3 subplots)
    ax_profiles = [fig.add_subplot(gs[0, i]) for i in range(3)]
    plot_concentration_profiles(results, ax_profiles)

    # Row 2 — Operating diagram (spans 2 cols) + legend/info panel
    ax_op = fig.add_subplot(gs[1, 0:2])
    plot_operating_diagram(results, ax_op)

    # Row 2, col 3 — Equations summary panel
    ax_eq = fig.add_subplot(gs[1, 2])
    ax_eq.set_facecolor(PANEL_COLOR)
    ax_eq.axis("off")
    eq_text = (
        "KEY EQUATIONS\n"
        "-------------------------\n\n"
        "Operating Lines (ODEs):\n"
        r"$\frac{dC_b}{dz}=\frac{1}{Q_b(z)}\left[-K_o\frac{A}{L}(C_b-C_d)"
        r"+\frac{Q_{uf}}{L}(1-S)C_b\right]$" + "\n\n"
        r"$\frac{dC_d}{dz}=\frac{1}{Q_d(z)}\left[-K_o\frac{A}{L}(C_b-C_d)"
        r"-\frac{Q_{uf}}{L}(SC_b-C_d)\right]$" + "\n\n"
        "Resistance-in-Series:\n"
        r"$\frac{1}{K_o}=\frac{1}{k_b}+\frac{\delta}{D_m}+\frac{1}{k_d}$"
        + "\n\n"
        "Stokes-Einstein:\n"
        r"$D = \frac{k_B T}{6\pi\mu\, r_s}$" + "\n\n"
        "Clearance:\n"
        r"$K = \frac{Q_{b,in}C_{b,in} - Q_{b,out}C_{b,out}}{C_{b,in}}$"
    )
    ax_eq.text(
        0.05, 0.95, eq_text, transform=ax_eq.transAxes,
        fontsize=10, color=TEXT_COLOR, verticalalignment="top",
        fontfamily="monospace", linespacing=1.6,
        bbox=dict(boxstyle="round,pad=0.5", facecolor=BG_COLOR,
                  edgecolor=ACCENT, alpha=0.9),
    )

    # Row 3 — Resistance pie charts
    ax_pies = [fig.add_subplot(gs[2, i]) for i in range(3)]
    plot_resistance_pies(qb, qd, ax_pies)

    return fig


# ══════════════════════════════════════════════════════════════════════════
#  SEPARATE EQUILIBRIUM / OPERATING-LINE FIGURE  (for close-up view)
# ══════════════════════════════════════════════════════════════════════════
def build_equilibrium_figure(results, qb, qd, quf):
    """
    Dedicated full-size operating diagram — the dialysis equivalent of
    a McCabe-Thiele y–x plot.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), facecolor=BG_COLOR)
    fig.suptitle(
        "Operating Lines vs Equilibrium Curve   |   "
        f"$Q_b$={qb:.0f}  $Q_d$={qd:.0f}  $Q_{{uf}}$={quf:.0f}  mL/min",
        fontsize=14, fontweight="bold", color=TEXT_COLOR, y=0.98,
    )

    # ── Left panel: Cb vs Cd (operating diagram) ──
    ax = axes[0]
    plot_operating_diagram(results, ax)

    # ── Right panel: Urea zoom-in with annotated driving force profile ──
    ax2 = axes[1]
    ax2.set_facecolor(PANEL_COLOR)

    if "Urea" in results:
        res = results["Urea"]
        z_cm = res["z"] * 100
        cb = res["Cb_profile"]
        cd = res["Cd_profile"]
        delta_c = cb - cd

        # Plot driving force along fiber length
        ax2.fill_between(z_cm, delta_c, alpha=0.25, color=COLORS["Urea"]["blood"])
        ax2.plot(z_cm, delta_c, lw=2.5, color=COLORS["Urea"]["blood"], label="$\\Delta C = C_b - C_d$")
        ax2.plot(z_cm, cb, lw=2, color=COLORS["Urea"]["blood"], ls="--", alpha=0.5, label="$C_b(z)$")
        ax2.plot(z_cm, cd, lw=2, color=COLORS["Urea"]["dial"], ls="--", alpha=0.5, label="$C_d(z)$")

        # Annotate min/max driving force
        i_max = np.argmax(delta_c)
        i_min = np.argmin(delta_c)
        ax2.annotate(f"dC_max = {delta_c[i_max]:.3f}",
                     xy=(z_cm[i_max], delta_c[i_max]),
                     xytext=(z_cm[i_max] + 2, delta_c[i_max] + 0.05),
                     fontsize=9, color=TEXT_COLOR, fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.5))
        ax2.annotate(f"dC_min = {delta_c[i_min]:.3f}",
                     xy=(z_cm[i_min], delta_c[i_min]),
                     xytext=(z_cm[i_min] - 5, delta_c[i_min] + 0.1),
                     fontsize=9, color=TEXT_COLOR, fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.5))

        # LMCD calculation
        dc1 = delta_c[0]
        dc2 = delta_c[-1]
        if abs(dc1 - dc2) > 1e-8:
            lmcd = (dc1 - dc2) / np.log(dc1 / dc2)
        else:
            lmcd = dc1
        ax2.axhline(lmcd, color=ACCENT, ls=":", lw=1.5, alpha=0.7)
        ax2.text(z_cm[-1] * 0.5, lmcd + 0.02,
                 f"LMCD = {lmcd:.4f} mg/mL", fontsize=10, color=ACCENT, fontweight="bold",
                 ha="center")

    ax2.set_xlabel("Fiber length z (cm)", fontsize=11, color=TEXT_COLOR)
    ax2.set_ylabel("Concentration (mg/mL)", fontsize=11, color=TEXT_COLOR)
    ax2.set_title("Urea — Driving Force Profile Along Fiber",
                  fontsize=13, fontweight="bold", color=TEXT_COLOR, pad=10)
    ax2.legend(fontsize=9, loc="best", framealpha=0.7,
               facecolor=PANEL_COLOR, edgecolor=GRID_COLOR, labelcolor=TEXT_COLOR)
    ax2.tick_params(colors=TEXT_COLOR, labelsize=9)
    ax2.grid(True, color=GRID_COLOR, linewidth=0.5, alpha=0.5)
    for spine in ax2.spines.values():
        spine.set_color(GRID_COLOR)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    return fig


# ══════════════════════════════════════════════════════════════════════════
#  INTERACTIVE MODE
# ══════════════════════════════════════════════════════════════════════════
def interactive_loop():
    """Let the presenter tweak Qb/Qd live and regenerate plots."""
    qb, qd, quf = 300.0, 500.0, 10.0
    plt.ion()

    print("\n+------------------------------------------------------+")
    print("|  INTERACTIVE DEMO MODE                              |")
    print("|  Type new values or press Enter to keep current.    |")
    print("|  Type 'q' to quit.                                  |")
    print("+------------------------------------------------------+\n")

    while True:
        print(f"  Current -> Qb={qb:.0f}, Qd={qd:.0f}, Quf={quf:.1f} mL/min")
        inp = input("  Enter Qb Qd [Quf]  (or 'q' to quit): ").strip()
        if inp.lower() in ("q", "quit", "exit"):
            break
        if inp:
            parts = inp.split()
            try:
                qb = float(parts[0])
                qd = float(parts[1])
                if len(parts) >= 3:
                    quf = float(parts[2])
            except (ValueError, IndexError):
                print("  [!] Enter numbers like: 250 600  or  250 600 15")
                continue

        plt.close("all")
        results = run_simulation(qb, qd, quf)
        print_results(results, qb, qd, quf)

        fig1 = build_demo_figure(results, qb, qd, quf)
        fig2 = build_equilibrium_figure(results, qb, qd, quf)
        plt.show(block=False)
        plt.pause(0.1)

    plt.close("all")
    print("\n  Demo ended. Good luck!\n")


# ══════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(description="Hemodialyzer Process Simulator Demo")
    parser.add_argument("--qb",  type=float, default=300.0, help="Blood flow rate (mL/min)")
    parser.add_argument("--qd",  type=float, default=500.0, help="Dialysate flow rate (mL/min)")
    parser.add_argument("--quf", type=float, default=10.0,  help="Ultrafiltration rate (mL/min)")
    parser.add_argument("--interactive", "-i", action="store_true",
                        help="Run in interactive mode (tweak params live)")
    parser.add_argument("--save", "-s", action="store_true",
                        help="Save plots to plots/ directory instead of showing")
    args = parser.parse_args()

    if args.interactive:
        interactive_loop()
        return

    # ── Single run ──
    results = run_simulation(args.qb, args.qd, args.quf)
    print_results(results, args.qb, args.qd, args.quf)

    # Build figures
    fig1 = build_demo_figure(results, args.qb, args.qd, args.quf)
    fig2 = build_equilibrium_figure(results, args.qb, args.qd, args.quf)

    if args.save:
        out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")
        os.makedirs(out_dir, exist_ok=True)
        fig1.savefig(os.path.join(out_dir, "demo_dashboard.png"), dpi=200, facecolor=BG_COLOR)
        fig2.savefig(os.path.join(out_dir, "operating_equilibrium.png"), dpi=200, facecolor=BG_COLOR)
        print(f"  [OK] Plots saved to {out_dir}/")
    else:
        plt.show()


if __name__ == "__main__":
    main()
