"""
Hemodialyzer Process Simulator — Interactive GUI
=================================================
Run:  streamlit run app.py
"""

import streamlit as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import sys, os

# ── path setup ───────────────────────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from parameters import SOLUTES, C_b_in, C_d_in
from solver import solve_solute_profile
from mass_transfer import calculate_overall_Ko

# ── theming constants ────────────────────────────────────────────────────
COLORS = {
    "Urea":    {"blood": "#E63946", "dial": "#457B9D", "op": "#E63946"},
    "VitB12":  {"blood": "#2A9D8F", "dial": "#264653", "op": "#2A9D8F"},
    "Albumin": {"blood": "#E9C46A", "dial": "#F4A261", "op": "#E9C46A"},
}
BG       = "#0F1117"
PANEL    = "#1A1B26"
GRID     = "#2A2B3D"
TXT      = "#E8E8F0"
ACCENT   = "#7AA2F7"
EQUIL    = "#FF6B6B"

plt.rcParams.update({
    "figure.facecolor": BG,
    "axes.facecolor":   PANEL,
    "axes.edgecolor":   GRID,
    "axes.labelcolor":  TXT,
    "text.color":       TXT,
    "xtick.color":      TXT,
    "ytick.color":      TXT,
    "grid.color":       GRID,
    "grid.alpha":       0.4,
})

# ═══════════════════════════════════════════════════════════════════
#  SOLVER
# ═══════════════════════════════════════════════════════════════════
@st.cache_data(show_spinner=False)
def run_sim(qb_ml, qd_ml, quf_ml):
    qb  = qb_ml  * 1e-6 / 60
    qd  = qd_ml  * 1e-6 / 60
    quf = quf_ml * 1e-6 / 60
    results = {}
    for s in SOLUTES:
        results[s] = solve_solute_profile(s, qb, qd, quf)
    return results

def resistance_breakdown(qb_ml, qd_ml):
    qb = qb_ml * 1e-6 / 60
    qd = qd_ml * 1e-6 / 60
    bd = {}
    for name, p in SOLUTES.items():
        Ko, Rb, Rm, Rd = calculate_overall_Ko(qb, qd, p["r_s"])
        Rt = Rb + Rm + Rd
        bd[name] = {"Ko": Ko, "Rb": 100*Rb/Rt, "Rm": 100*Rm/Rt, "Rd": 100*Rd/Rt}
    return bd

# ═══════════════════════════════════════════════════════════════════
#  PLOT BUILDERS
# ═══════════════════════════════════════════════════════════════════
def fig_profiles(results):
    """Row of 3 concentration-profile plots."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    for ax, (name, res) in zip(axes, results.items()):
        z = res["z"] * 100
        cb, cd = res["Cb_profile"], res["Cd_profile"]
        ax.fill_between(z, cb, alpha=0.15, color=COLORS[name]["blood"])
        ax.fill_between(z, cd, alpha=0.10, color=COLORS[name]["dial"])
        ax.plot(z, cb, lw=2.5, color=COLORS[name]["blood"], label="$C_b(z)$ Blood")
        ax.plot(z, cd, lw=2.5, color=COLORS[name]["dial"], ls="--", label="$C_d(z)$ Dialysate")
        ax.set_title(name, fontsize=13, fontweight="bold", pad=8)
        ax.set_xlabel("Fiber length z (cm)", fontsize=10)
        ax.set_ylabel("Conc. (mg/mL)", fontsize=10)
        ax.legend(fontsize=8, loc="best", framealpha=0.5, facecolor=PANEL,
                  edgecolor=GRID, labelcolor=TXT)
        ax.grid(True, linewidth=0.5)
        K = res["Clearance_ml_min"]
        ax.text(0.97, 0.95, f"K = {K:.1f} mL/min", transform=ax.transAxes,
                ha="right", va="top", fontsize=10, fontweight="bold", color=ACCENT,
                bbox=dict(boxstyle="round,pad=0.3", fc=BG, ec=ACCENT, alpha=0.85))
    fig.tight_layout()
    return fig


def fig_operating(results):
    """Cb-vs-Cd operating diagram with 45-degree equilibrium line."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))

    # ── left: Cb vs Cd ──
    ax = axes[0]
    max_c = 0
    handles = []
    for name, res in results.items():
        cb, cd = res["Cb_profile"], res["Cd_profile"]
        color = COLORS[name]["op"]
        if name == "Albumin" and res["Clearance_ml_min"] < 0.5:
            continue
        max_c = max(max_c, cb.max(), cd.max())
        ax.plot(cd, cb, lw=2.5, color=color, zorder=3)
        ax.scatter(cd[0],  cb[0],  s=80, color=color, zorder=5,
                   edgecolors="white", linewidths=1.2)
        ax.scatter(cd[-1], cb[-1], s=80, color=color, zorder=5,
                   edgecolors="white", linewidths=1.2, marker="s")
        mid = len(cb)//2
        df = cb[mid] - cd[mid]
        ax.annotate("", xy=(cb[mid], cb[mid]), xytext=(cd[mid], cb[mid]),
                    arrowprops=dict(arrowstyle="<->", color=color, lw=1.5, alpha=0.7))
        ax.text((cd[mid]+cb[mid])/2, cb[mid]+max_c*0.02,
                f"dC={df:.3f}", ha="center", fontsize=8, color=color, fontweight="bold")
        handles.append(Line2D([0],[0], color=color, lw=2.5, label=name))

    lim = max_c * 1.15
    ax.plot([0, lim], [0, lim], ls="--", lw=2, color=EQUIL, alpha=0.8, zorder=2)
    handles.append(Line2D([0],[0], color=EQUIL, lw=2, ls="--",
                          label="Equilibrium ($C_b=C_d$)"))

    if "Urea" in results:
        cb_u = results["Urea"]["Cb_profile"]
        cd_u = results["Urea"]["Cd_profile"]
        idx = np.argsort(cd_u)
        ax.fill_betweenx(cb_u[idx], cd_u[idx], cb_u[idx],
                         alpha=0.08, color=COLORS["Urea"]["op"])

    ax.set_xlim(left=-0.02)
    ax.set_ylim(bottom=-0.02)
    ax.set_xlabel("$C_d$ - Dialysate (mg/mL)", fontsize=11)
    ax.set_ylabel("$C_b$ - Blood (mg/mL)", fontsize=11)
    ax.set_title("Operating Diagram (Cb vs Cd)\nEquilibrium = 45 degree line, gap = driving force",
                 fontsize=12, fontweight="bold", pad=10)
    ax.legend(handles=handles, fontsize=9, loc="upper left", framealpha=0.6,
              facecolor=PANEL, edgecolor=GRID, labelcolor=TXT)
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linewidth=0.5)

    # ── right: driving-force profile for Urea ──
    ax2 = axes[1]
    if "Urea" in results:
        res = results["Urea"]
        z = res["z"] * 100
        cb, cd = res["Cb_profile"], res["Cd_profile"]
        dc = cb - cd
        ax2.fill_between(z, dc, alpha=0.25, color=COLORS["Urea"]["blood"])
        ax2.plot(z, dc, lw=2.5, color=COLORS["Urea"]["blood"],
                 label=r"$\Delta C = C_b - C_d$")
        ax2.plot(z, cb, lw=2, color=COLORS["Urea"]["blood"], ls="--", alpha=0.5,
                 label="$C_b(z)$")
        ax2.plot(z, cd, lw=2, color=COLORS["Urea"]["dial"], ls="--", alpha=0.5,
                 label="$C_d(z)$")
        # LMCD
        dc1, dc2 = dc[0], dc[-1]
        lmcd = (dc1 - dc2) / np.log(dc1/dc2) if abs(dc1 - dc2) > 1e-8 else dc1
        ax2.axhline(lmcd, color=ACCENT, ls=":", lw=1.5, alpha=0.7)
        ax2.text(z[-1]*0.5, lmcd + 0.025, f"LMCD = {lmcd:.4f} mg/mL",
                 fontsize=10, color=ACCENT, fontweight="bold", ha="center")
        # annotate extremes
        i_mx, i_mn = np.argmax(dc), np.argmin(dc)
        ax2.annotate(f"dC_max = {dc[i_mx]:.3f}", xy=(z[i_mx], dc[i_mx]),
                     xytext=(z[i_mx]+2, dc[i_mx]+0.05), fontsize=9,
                     fontweight="bold", color=TXT,
                     arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.5))
        ax2.annotate(f"dC_min = {dc[i_mn]:.3f}", xy=(z[i_mn], dc[i_mn]),
                     xytext=(z[i_mn]-6, dc[i_mn]+0.12), fontsize=9,
                     fontweight="bold", color=TXT,
                     arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.5))
    ax2.set_xlabel("Fiber length z (cm)", fontsize=11)
    ax2.set_ylabel("Concentration (mg/mL)", fontsize=11)
    ax2.set_title("Urea - Driving Force Along Fiber", fontsize=12,
                  fontweight="bold", pad=10)
    ax2.legend(fontsize=9, loc="best", framealpha=0.6, facecolor=PANEL,
               edgecolor=GRID, labelcolor=TXT)
    ax2.grid(True, linewidth=0.5)
    fig.tight_layout()
    return fig


def fig_pies(qb, qd):
    """Resistance breakdown pie charts."""
    bd = resistance_breakdown(qb, qd)
    pie_c = ["#E63946", "#457B9D", "#2A9D8F"]
    labels = ["Blood film", "Membrane", "Dialysate film"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (name, info) in zip(axes, bd.items()):
        sizes = [info["Rb"], info["Rm"], info["Rd"]]
        wedges, texts, autos = ax.pie(
            sizes, labels=labels, autopct="%1.1f%%", startangle=140,
            colors=pie_c, textprops={"fontsize": 8, "color": TXT},
            wedgeprops={"edgecolor": BG, "linewidth": 1.5}, pctdistance=0.6)
        for t in autos:
            t.set_fontsize(8); t.set_fontweight("bold")
        ax.set_title(f"{name}\n$K_o$ = {info['Ko']:.2e} m/s",
                     fontsize=11, fontweight="bold", pad=6)
    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════
#  STREAMLIT APP
# ═══════════════════════════════════════════════════════════════════
st.set_page_config(
    page_title="Hemodialyzer Simulator",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── custom CSS ───────────────────────────────────────────────────
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');

    html, body, [class*="css"] { font-family: 'Inter', sans-serif; }

    /* ─── global dark theme ─── */
    .stApp { background: #0F1117; color: #E8E8F0 !important; }

    /* force ALL text to be visible */
    p, span, label, li, td, th, div, a,
    .stMarkdown, .stMarkdown p, .stMarkdown span,
    [data-testid="stMarkdownContainer"],
    [data-testid="stMarkdownContainer"] p,
    [data-testid="stMarkdownContainer"] span,
    [data-testid="stMarkdownContainer"] li {
        color: #E8E8F0 !important;
    }

    /* sidebar */
    section[data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1A1B26 0%, #0F1117 100%);
        border-right: 1px solid #2A2B3D;
    }
    section[data-testid="stSidebar"] h2,
    section[data-testid="stSidebar"] h3,
    section[data-testid="stSidebar"] .stMarkdown h2,
    section[data-testid="stSidebar"] .stMarkdown h3 {
        color: #7AA2F7 !important;
    }
    section[data-testid="stSidebar"] p,
    section[data-testid="stSidebar"] span,
    section[data-testid="stSidebar"] label,
    section[data-testid="stSidebar"] td,
    section[data-testid="stSidebar"] th,
    section[data-testid="stSidebar"] li {
        color: #E8E8F0 !important;
    }

    /* slider labels & values */
    .stSlider label, .stSlider span,
    [data-testid="stSlider"] label,
    [data-testid="stSlider"] span,
    .stSlider [data-testid="stThumbValue"],
    [data-testid="stThumbValue"] {
        color: #E8E8F0 !important;
    }

    /* tab labels */
    .stTabs [data-baseweb="tab-list"] button {
        color: #9CA3AF !important;
    }
    .stTabs [data-baseweb="tab-list"] button[aria-selected="true"] {
        color: #7AA2F7 !important;
        border-bottom-color: #7AA2F7 !important;
    }

    /* blockquotes */
    blockquote, blockquote p, blockquote span,
    .stMarkdown blockquote, .stMarkdown blockquote p {
        color: #C0CAF5 !important;
    }

    /* metric cards */
    div[data-testid="stMetric"] {
        background: linear-gradient(135deg, #1A1B26, #252736);
        border: 1px solid #2A2B3D;
        border-radius: 12px;
        padding: 16px 20px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.3);
    }
    div[data-testid="stMetric"] label {
        color: #9CA3AF !important;
        font-size: 0.85rem !important;
    }
    div[data-testid="stMetric"] [data-testid="stMetricValue"] {
        color: #7AA2F7 !important;
        font-weight: 700 !important;
        font-size: 1.6rem !important;
    }

    /* dividers */
    hr { border-color: #2A2B3D !important; }

    /* header text */
    h1 { color: #E8E8F0 !important; }
    h2, h3 { color: #C0CAF5 !important; }
    h4 { color: #7AA2F7 !important; }

    /* markdown tables in sidebar etc */
    .stMarkdown table { color: #E8E8F0 !important; }
    .stMarkdown table th { color: #7AA2F7 !important; background: #252736 !important; }
    .stMarkdown table td { color: #E8E8F0 !important; border-color: #2A2B3D !important; }

    /* info boxes */
    .equation-box {
        background: linear-gradient(135deg, #1A1B26, #1E2030);
        border: 1px solid #7AA2F7;
        border-radius: 10px;
        padding: 20px;
        margin: 10px 0;
        color: #E8E8F0 !important;
        font-family: 'Inter', sans-serif;
    }
    .equation-box h4 { color: #7AA2F7 !important; margin-top: 0; }
    .equation-box p { color: #E8E8F0 !important; }

    /* styled results table */
    .styled-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.95rem;
        color: #E8E8F0 !important;
    }
    .styled-table th {
        background: #252736;
        color: #7AA2F7 !important;
        padding: 12px 16px;
        text-align: left;
        border-bottom: 2px solid #7AA2F7;
        font-weight: 600;
    }
    .styled-table td {
        padding: 10px 16px;
        border-bottom: 1px solid #2A2B3D;
        color: #E8E8F0 !important;
    }
    .styled-table tr:hover td { background: #1E2030; }

    /* tooltips / help icons */
    [data-testid="stTooltipIcon"] { color: #7AA2F7 !important; }

    /* hide streamlit branding but keep sidebar toggle */
    #MainMenu { visibility: hidden; }
    footer { visibility: hidden; }
    header [data-testid="stHeader"] { background: #0F1117 !important; }
</style>
""", unsafe_allow_html=True)

# ── sidebar controls ─────────────────────────────────────────────
with st.sidebar:
    st.markdown("## :control_knobs: Operating Conditions")
    st.markdown("---")
    qb = st.slider("Blood Flow Rate  Qb (mL/min)", 100, 500, 300, step=10,
                    help="Typical clinical range: 200-450 mL/min")
    qd = st.slider("Dialysate Flow Rate  Qd (mL/min)", 200, 800, 500, step=10,
                    help="Typical clinical range: 300-800 mL/min")
    quf = st.slider("Ultrafiltration  Quf (mL/min)", 0, 30, 10, step=1,
                     help="Net fluid removal rate")

    st.markdown("---")
    st.markdown("## :memo: Membrane Parameters")
    st.markdown(f"""
    | Parameter | Value |
    |---|---|
    | Fiber ID | 200 um |
    | Membrane thickness | 40 um |
    | Fiber length | 25 cm |
    | Number of fibers | 10,000 |
    | Porosity | 0.75 |
    | Tortuosity | 2.5 |
    | Pore radius | 3.0 nm |
    """)

    st.markdown("---")
    st.markdown(
        "<div style='text-align:center; color:#555; font-size:0.8rem;'>"
        "Separation Processes Course Project<br>Apr 2026</div>",
        unsafe_allow_html=True,
    )

# ── main content ─────────────────────────────────────────────────
st.markdown(
    "<h1 style='text-align:center; margin-bottom:0;'>"
    "Hemodialyzer Process Simulator</h1>"
    "<p style='text-align:center; color:#9CA3AF; font-size:1.1rem; margin-top:4px;'>"
    "Counter-Current Hollow-Fiber Membrane Mass Transfer</p>",
    unsafe_allow_html=True
)

# run solver
with st.spinner("Solving ODEs..."):
    results = run_sim(qb, qd, quf)

# ── metrics row ──────────────────────────────────────────────────
st.markdown("---")
cols = st.columns(5)
cols[0].metric("Blood Flow (Qb)", f"{qb} mL/min")
cols[1].metric("Dialysate (Qd)", f"{qd} mL/min")

solute_labels = {"Urea": "Urea K", "VitB12": "VitB12 K", "Albumin": "Albumin K"}
for i, (name, res) in enumerate(results.items()):
    K = res["Clearance_ml_min"]
    cols[i+2].metric(solute_labels[name], f"{K:.1f} mL/min")

# ── tabs ─────────────────────────────────────────────────────────
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Concentration Profiles",
    "Operating Diagram & Equilibrium",
    "Resistance Breakdown",
    "Equations & Theory",
    "ML Surrogate Model",
])

with tab1:
    st.markdown("### Spatial Concentration Profiles Along Fiber Length")
    st.markdown(
        "> Blood enters at z = 0, dialysate enters at z = L (counter-current). "
        "Urea transfers rapidly; Albumin is fully rejected."
    )
    st.pyplot(fig_profiles(results), use_container_width=True)

with tab2:
    st.markdown("### Operating Lines vs Equilibrium Curve")
    st.markdown(
        "> **Left:** Cb vs Cd plot -- the dialysis analog of a McCabe-Thiele diagram. "
        "The 45-degree dashed line is the equilibrium (Cb = Cd, no driving force). "
        "The gap between the operating curve and this line is the local driving force.  \n"
        "> **Right:** Driving force profile for Urea along fiber length, with LMCD annotated."
    )
    st.pyplot(fig_operating(results), use_container_width=True)

with tab3:
    st.markdown("### Mass Transfer Resistance Breakdown")
    st.markdown(
        "> For **Urea** (small molecule), blood-side film + membrane share resistance.  \n"
        "> For **VitB12** (middle molecule), membrane dominates due to steric hindrance.  \n"
        "> For **Albumin** (large protein), membrane resistance is 100% -- fully rejected."
    )
    st.pyplot(fig_pies(qb, qd), use_container_width=True)

    # ── results table ──
    st.markdown("### Detailed Results")
    bd = resistance_breakdown(qb, qd)
    table_html = """<table class="styled-table"><thead><tr>
        <th>Solute</th><th>MW (Da)</th><th>Sieving Coeff</th>
        <th>Clearance (mL/min)</th><th>Cb_out (mg/mL)</th>
        <th>Ko (m/s)</th><th>Blood %</th><th>Membrane %</th><th>Dialysate %</th>
    </tr></thead><tbody>"""
    for name, res in results.items():
        info = bd[name]
        table_html += (
            f"<tr><td><strong>{name}</strong></td>"
            f"<td>{SOLUTES[name]['Mw']:.0f}</td>"
            f"<td>{res['Sieving_Coefficient']:.4f}</td>"
            f"<td>{res['Clearance_ml_min']:.1f}</td>"
            f"<td>{res['Cb_out']:.4f}</td>"
            f"<td>{info['Ko']:.2e}</td>"
            f"<td>{info['Rb']:.1f}%</td>"
            f"<td>{info['Rm']:.1f}%</td>"
            f"<td>{info['Rd']:.1f}%</td></tr>"
        )
    table_html += "</tbody></table>"
    st.markdown(table_html, unsafe_allow_html=True)

with tab4:
    st.markdown("### Key Equations")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("""
<div class="equation-box">
<h4>Operating Lines (Coupled ODEs)</h4>
<p>Blood side:</p>

$$\\frac{dC_b}{dz} = \\frac{1}{Q_b(z)}\\left[-K_o\\frac{A}{L}(C_b - C_d) + \\frac{Q_{uf}}{L}(1-S)C_b\\right]$$

<p>Dialysate side:</p>

$$\\frac{dC_d}{dz} = \\frac{1}{Q_d(z)}\\left[-K_o\\frac{A}{L}(C_b - C_d) - \\frac{Q_{uf}}{L}(SC_b - C_d)\\right]$$

<p style="color:#9CA3AF; font-size:0.85rem;">
Solved as a Boundary Value Problem: C<sub>b</sub>(0) = C<sub>b,in</sub> and C<sub>d</sub>(L) = 0
</p>
</div>

<div class="equation-box">
<h4>Stokes-Einstein Diffusivity</h4>

$$D = \\frac{k_B T}{6\\pi\\mu\\, r_s}$$

</div>

<div class="equation-box">
<h4>Membrane Diffusivity (with steric hindrance)</h4>

$$D_m = D_{water} \\cdot \\frac{\\varepsilon}{\\tau} \\cdot \\left(1 - \\frac{r_s}{r_{pore}}\\right)^2$$

</div>
""", unsafe_allow_html=True)

    with c2:
        st.markdown("""
<div class="equation-box">
<h4>Resistance-in-Series Model</h4>

$$\\frac{1}{K_o} = \\frac{1}{k_b} + \\frac{\\delta}{D_m} + \\frac{1}{k_d}$$

<p style="color:#9CA3AF; font-size:0.85rem;">
k<sub>b</sub> from Leveque (Sh = 1.62 Gz<sup>1/3</sup>),
k<sub>d</sub> from shell-side Sherwood correlation
</p>
</div>

<div class="equation-box">
<h4>Sieving Coefficient (Ferry equation)</h4>

$$S = (1-\\lambda)^2 (2 - (1-\\lambda)^2), \\quad \\lambda = r_s / r_{pore}$$

</div>

<div class="equation-box">
<h4>Clearance</h4>

$$K = \\frac{Q_{b,in} C_{b,in} - Q_{b,out} C_{b,out}}{C_{b,in}}$$

</div>

<div class="equation-box">
<h4>Log Mean Concentration Difference</h4>

$$LMCD = \\frac{\\Delta C_1 - \\Delta C_2}{\\ln(\\Delta C_1 / \\Delta C_2)}$$

</div>
""", unsafe_allow_html=True)

with tab5:
    st.markdown("### Digital Twin — ML Surrogate Model")
    st.markdown(
        "> A **Random Forest** is trained on 1,000 ODE solver runs to predict clearance instantly. "
        "This acts as a **digital twin** — a fast surrogate that replaces the expensive BVP solver."
    )

    # ── train / load model ─────────────────────────────────────────
    import pandas as pd
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import (mean_absolute_error,
                                 root_mean_squared_error,
                                 mean_absolute_percentage_error)

    @st.cache_resource(show_spinner="Training ML surrogate model...")
    def get_ml_results():
        csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "data", "dialyzer_dataset.csv")
        if not os.path.exists(csv_path):
            st.error("Dataset not found. Run `python src/ml_model.py` first.")
            return None

        df = pd.read_csv(csv_path)
        X = df[["Qb_ml_min", "Qd_ml_min", "Quf_ml_min"]]
        y = df[["Clearance_Urea", "Clearance_B12"]]
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42)

        model = RandomForestRegressor(n_estimators=100, random_state=42)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)

        return {
            "model": model,
            "X_test": X_test, "y_test": y_test, "y_pred": y_pred,
            "features": X.columns.tolist(),
        }

    ml = get_ml_results()

    if ml is not None:
        y_test = ml["y_test"]
        y_pred = ml["y_pred"]

        # ── compute metrics ──
        solutes_ml = ["Urea", "VitB12"]
        metrics = {}
        for i, sol in enumerate(solutes_ml):
            yt = y_test.iloc[:, i].values
            yp = y_pred[:, i]
            metrics[sol] = {
                "MAE":  mean_absolute_error(yt, yp),
                "RMSE": root_mean_squared_error(yt, yp),
                "MAPE": mean_absolute_percentage_error(yt, yp) * 100,
                "Max Error": np.max(np.abs(yt - yp)),
                "yt": yt, "yp": yp,
            }

        # ── metric cards ──
        st.markdown("#### Model Performance Metrics")
        mc1, mc2 = st.columns(2)
        for col, sol in zip([mc1, mc2], solutes_ml):
            m = metrics[sol]
            with col:
                st.markdown(f"**{sol} Clearance**")
                r1, r2, r3, r4 = st.columns(4)
                r1.metric("MAE", f"{m['MAE']:.2f} mL/min")
                r2.metric("RMSE", f"{m['RMSE']:.2f} mL/min")
                r3.metric("MAPE", f"{m['MAPE']:.2f}%")
                r4.metric("Max Err", f"{m['Max Error']:.2f} mL/min")

        st.markdown("---")

        # ── 1. Actual vs Predicted ──
        st.markdown("#### Actual vs Predicted Clearance")
        fig_avp, axes_avp = plt.subplots(1, 2, figsize=(14, 5.5))
        for ax, sol in zip(axes_avp, solutes_ml):
            m = metrics[sol]
            yt, yp = m["yt"], m["yp"]
            color = COLORS.get(sol, COLORS["Urea"])["blood"]
            ax.scatter(yt, yp, alpha=0.5, s=25, color=color, edgecolors="white",
                       linewidths=0.3, zorder=3)
            mn, mx = min(yt.min(), yp.min()), max(yt.max(), yp.max())
            pad = (mx - mn) * 0.05
            ax.plot([mn-pad, mx+pad], [mn-pad, mx+pad], "--", color=EQUIL,
                    lw=2, alpha=0.8, label="Perfect prediction")
            ax.set_xlabel("ODE Solver (Actual) [mL/min]", fontsize=10)
            ax.set_ylabel("Random Forest (Predicted) [mL/min]", fontsize=10)
            ax.set_title(f"{sol} — MAE={m['MAE']:.2f}, RMSE={m['RMSE']:.2f}",
                         fontsize=12, fontweight="bold")
            ax.legend(fontsize=9, loc="upper left", framealpha=0.6,
                      facecolor=PANEL, edgecolor=GRID, labelcolor=TXT)
            ax.grid(True, linewidth=0.5)
        fig_avp.tight_layout()
        st.pyplot(fig_avp, use_container_width=True)

        # ── 2. Residual Distribution ──
        st.markdown("#### Prediction Error Distribution")
        fig_res, axes_res = plt.subplots(1, 2, figsize=(14, 4.5))
        for ax, sol in zip(axes_res, solutes_ml):
            m = metrics[sol]
            errors = m["yt"] - m["yp"]
            color = COLORS.get(sol, COLORS["Urea"])["blood"]
            ax.hist(errors, bins=30, color=color, alpha=0.7,
                    edgecolor=BG, linewidth=0.8)
            ax.axvline(0, color=EQUIL, ls="--", lw=2, alpha=0.8)
            ax.axvline(errors.mean(), color=ACCENT, ls=":", lw=2, alpha=0.8)
            ax.text(0.97, 0.95,
                    f"Mean err = {errors.mean():.3f}\nStd = {errors.std():.3f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=9,
                    fontweight="bold", color=TXT,
                    bbox=dict(boxstyle="round,pad=0.3", fc=BG, ec=ACCENT, alpha=0.85))
            ax.set_xlabel("Error: Actual - Predicted (mL/min)", fontsize=10)
            ax.set_ylabel("Count", fontsize=10)
            ax.set_title(f"{sol} — Residual Distribution",
                         fontsize=12, fontweight="bold")
            ax.grid(True, linewidth=0.5)
        fig_res.tight_layout()
        st.pyplot(fig_res, use_container_width=True)

        # ── 3. Feature Importance ──
        st.markdown("#### Feature Importance")
        importances = ml["model"].feature_importances_
        feature_labels = ["Blood Flow (Qb)", "Dialysate Flow (Qd)",
                          "Ultrafiltration (Quf)"]

        # get per-target importance if available
        fig_imp, ax_imp = plt.subplots(figsize=(10, 4))
        x_pos = np.arange(len(feature_labels))
        bar_colors = ["#E63946", "#457B9D", "#2A9D8F"]

        # average importance across both targets
        per_target_imp = np.array([tree.feature_importances_
                                   for tree in ml["model"].estimators_])
        mean_imp = per_target_imp.mean(axis=0)
        std_imp = per_target_imp.std(axis=0)

        bars = ax_imp.bar(x_pos, mean_imp, yerr=std_imp, capsize=6,
                          color=bar_colors, edgecolor=BG, linewidth=1.5,
                          alpha=0.85, zorder=3)
        for bar, val in zip(bars, mean_imp):
            ax_imp.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                        f"{val:.3f}", ha="center", va="bottom", fontsize=11,
                        fontweight="bold", color=TXT)
        ax_imp.set_xticks(x_pos)
        ax_imp.set_xticklabels(feature_labels, fontsize=11)
        ax_imp.set_ylabel("Importance (mean +/- std across trees)", fontsize=10)
        ax_imp.set_title("Random Forest Feature Importance",
                         fontsize=13, fontweight="bold")
        ax_imp.set_ylim(0, max(mean_imp) * 1.25)
        ax_imp.grid(True, axis="y", linewidth=0.5)
        fig_imp.tight_layout()
        st.pyplot(fig_imp, use_container_width=True)

        # ── explanation ──
        st.markdown("""
<div class="equation-box">
<h4>How the Surrogate Model Works</h4>
<p>1. The rigorous ODE solver was run <strong>1,000 times</strong> with randomly sampled operating
conditions (Qb: 200-500, Qd: 300-800, Quf: 0-20 mL/min).</p>
<p>2. A <strong>Random Forest Regressor</strong> (100 trees) was trained on 80% of the data to predict
Urea and VitB12 clearance from the three operating parameters.</p>
<p>3. The model is evaluated on the held-out 20% test set using:</p>
<table class="styled-table" style="margin:10px 0;">
<thead><tr><th>Metric</th><th>What It Measures</th><th>Why Not R&sup2;?</th></tr></thead>
<tbody>
<tr><td><strong>MAE</strong></td><td>Average absolute error in mL/min</td><td>Directly interpretable in physical units</td></tr>
<tr><td><strong>RMSE</strong></td><td>Root mean squared error &mdash; penalizes large errors</td><td>Sensitive to outliers, more conservative</td></tr>
<tr><td><strong>MAPE</strong></td><td>Percentage error relative to actual value</td><td>Scale-independent, shows relative accuracy</td></tr>
<tr><td><strong>Max Error</strong></td><td>Worst-case single prediction error</td><td>Important for safety-critical applications</td></tr>
</tbody>
</table>
<p style="color:#9CA3AF; font-size:0.85rem;">
The ODE solver takes ~0.5s per run. The ML surrogate predicts in microseconds &mdash;
enabling real-time optimization and control applications.
</p>
</div>
""", unsafe_allow_html=True)

