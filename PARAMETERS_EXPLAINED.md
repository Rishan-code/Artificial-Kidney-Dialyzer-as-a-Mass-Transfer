# Understanding the Parameters

Every parameter in this simulator represents a real, measurable physical quantity. This document explains what each one is, where it comes from, and what happens when you change it.

---

## 1. Solute Properties

These define the three molecules we track as they transfer from blood to dialysate.

### Why these three?

They represent the full spectrum of molecules in blood — from tiny waste (urea) to essential proteins (albumin). The dialyzer must remove the small ones and reject the large ones.

| Parameter | Urea | Vitamin B12 | Albumin | Unit |
|---|---|---|---|---|
| **Molecular Weight (Mw)** | 60.06 | 1,355 | 66,000 | g/mol (Da) |
| **Stokes Radius (r_s)** | 2.6 × 10⁻¹⁰ | 8.5 × 10⁻¹⁰ | 35.8 × 10⁻¹⁰ | m |
| **Inlet Blood Conc (C_b_in)** | 1.5 | 0.02 | 40.0 | mg/mL |
| **Clinical role** | Primary uremic toxin | Middle molecule marker | Essential blood protein |  |

### Molecular Weight (Mw)

The mass of one mole of molecules. Larger Mw generally means larger molecular size.

- **Urea (60 Da):** The main waste product from protein metabolism. Healthy kidneys excrete ~30g/day. In kidney failure, it accumulates to toxic levels.
- **VitB12 (1,355 Da):** Used as a representative "middle molecule." In practice, beta-2 microglobulin (11,800 Da) is more clinically relevant, but VitB12 is a standard test marker.
- **Albumin (66,000 Da):** The most abundant blood protein. Carries hormones, drugs, and maintains oncotic pressure. Losing albumin = bad outcome.

### Stokes Radius (r_s)

The effective hydrodynamic radius of the molecule — how large it appears to the surrounding fluid as it diffuses. Measured experimentally through diffusion experiments.

**Where it matters:**
1. **Diffusivity:** D = kT / (6πμ·r_s). Larger r_s → slower diffusion.
2. **Steric hindrance:** λ = r_s / r_pore. If λ ≥ 1, the molecule is blocked entirely.
3. **Sieving coefficient:** S = (1-λ)²(2-(1-λ)²). Determines convective transport.

| Solute | r_s (nm) | λ = r_s / r_pore | Sieving (S) | Passes through? |
|---|---|---|---|---|
| Urea | 0.26 | 0.087 | 0.97 | Yes, freely |
| VitB12 | 0.85 | 0.283 | 0.76 | Yes, with friction |
| Albumin | 3.58 | 1.193 | 0.00 | No — bigger than pore |

### Inlet Concentrations (C_b_in)

The concentration in the patient's blood before dialysis begins.

- **Urea = 1.5 mg/mL:** Corresponds to BUN (Blood Urea Nitrogen) of ~70 mg/dL, typical for a pre-dialysis patient. Normal range: 7-20 mg/dL. This patient is ~3.5x above normal.
- **VitB12 = 0.02 mg/mL:** Representative concentration for simulation purposes.
- **Albumin = 40 mg/mL:** Equivalent to 4.0 g/dL, right in the normal range (3.5-5.0 g/dL).
- **Dialysate inlet (C_d_in) = 0:** Fresh dialysate contains no waste. This maximizes the concentration gradient for diffusion.

---

## 2. Membrane Structural Parameters

These define the physical construction of the hollow fiber membrane.

```
    Cross-section of ONE fiber:
    
         <---- d_i = 200 μm ---->
         <-- delta = 40 μm -->
    ┌────┬────────────────────┬────┐
    │wall│     blood lumen     │wall│
    │    │    (blood flows     │    │
    │    │     through here)   │    │  
    └────┴────────────────────┴────┘
    
    10,000 of these fibers packed in parallel
    Each fiber is L = 25 cm long
```

### Inner Diameter (d_i = 200 μm)

The internal diameter of each hollow fiber — where blood flows.

- **Realistic range:** 180-220 μm for commercial dialyzers
- **Effect:** Smaller d_i → higher blood velocity → thinner boundary layer → better kb. But too small → excessive pressure drop and risk of clotting.
- **200 μm** is the industry standard balance.

### Membrane Thickness (δ = 40 μm)

The wall thickness of each fiber — the barrier between blood and dialysate.

- **Realistic range:** 30-50 μm
- **Effect:** Thinner membrane → less resistance (δ/Dm is smaller) → better Ko. But too thin → mechanical fragility, risk of rupture.
- **Appears in:** 1/Ko = 1/kb + **δ/Dm** + 1/kd

### Fiber Length (L = 25 cm)

How long each fiber is, from blood inlet to blood outlet.

- **Realistic range:** 20-30 cm
- **Effect:** Longer fiber → more membrane area → more total transfer. But also → higher pressure drop and potential for boundary layer development (reduces kb at the outlet end).

### Number of Fibers (n_fibers = 10,000)

Total number of parallel fibers packed into the dialyzer housing.

- **Realistic range:** 8,000-16,000 depending on dialyzer model
- **Effect:** More fibers → more total membrane area (A = n × π × d_i × L) → more transfer capacity. The total area here is: 10,000 × π × 200×10⁻⁶ × 0.25 = **1.57 m²** (typical for a mid-size dialyzer).
- Each fiber carries: 300/10,000 = 0.03 mL/min of blood

### Porosity (ε = 0.75)

The void fraction of the membrane — what percentage is open pore space vs solid polymer.

- **Realistic range:** 0.60-0.80 for polysulfone membranes
- **Effect:** Higher porosity → more pathways for molecules → higher Dm → lower membrane resistance.
- **Appears in:** D_m = D_water × **(ε/τ)** × (1 - r_s/r_pore)²

### Tortuosity (τ = 2.5)

How winding the pore channels are. A straight-through pore has τ = 1. Real membrane pores zigzag through the wall.

- **Realistic range:** 2.0-3.0
- **Effect:** Higher tortuosity → longer effective path → lower Dm → higher membrane resistance.
- **Physical meaning:** τ = 2.5 means the actual path a molecule travels through the membrane is 2.5× longer than the membrane thickness.

### Pore Radius (r_pore = 3.0 nm)

The average radius of the pores in the membrane. This is the **most critical design parameter** — it determines what passes through and what doesn't.

- **Realistic range:** 2-5 nm for "high-flux" dialyzers; <2 nm for "low-flux"
- **Effect:** Larger pores → more molecules pass through (including potentially harmful protein loss). Smaller pores → better selectivity but lower clearance for middle molecules.
- **3 nm** provides: Urea passes freely, VitB12 passes with hindrance, Albumin blocked.

---

## 3. Patient & Fluid Parameters

### Temperature (T = 37°C = 310.15 K)

Body temperature. Affects:
- **Diffusivity** via Stokes-Einstein: D ∝ T. Higher T → faster diffusion (molecules have more kinetic energy).
- **Viscosity** of water: lower at higher T → also increases D.

In practice, dialysate is warmed to 37°C to match body temperature and prevent hypothermia.

### Hematocrit (Hct = 0.40)

The fraction of blood volume occupied by red blood cells (40%).

- **Realistic range:** 0.33-0.45 (lower in dialysis patients due to anemia)
- **Effect:** Used to calculate blood viscosity: **μ_blood = μ_water × (1 + 2.5 × Hct)**
  - At Hct = 0.40: μ_blood = 6.9×10⁻⁴ × 2.0 = 1.38×10⁻³ Pa·s
  - Blood is ~2× more viscous than water due to suspended red blood cells
- Higher Hct → higher viscosity → lower Reynolds number → thicker boundary layer → lower kb

### Water Viscosity (μ_water = 6.9 × 10⁻⁴ Pa·s)

Dynamic viscosity of water at 37°C. This is a well-established physical constant.

- Used in Stokes-Einstein: D = kBT / (6π × **μ_water** × r_s)
- Used to estimate blood viscosity via Einstein's viscosity equation

### Boltzmann Constant (kB = 1.38 × 10⁻²³ J/K)

Fundamental physical constant relating energy to temperature at the molecular level.

- Used only in Stokes-Einstein diffusivity calculation.

---

## 4. Operating Parameters (The Knobs You Control)

These are the three parameters you adjust via the sliders in the simulator.

### Blood Flow Rate (Qb = 300 mL/min)

The rate at which the blood pump pushes blood through the fibers.

- **Clinical range:** 200-450 mL/min
- **Typical prescription:** 300-400 mL/min
- **What it affects:**
  - Blood velocity in each fiber: v = Qb / (n × π × (d_i/2)²)
  - Reynolds number → Sherwood number → **blood-side mass transfer coefficient (kb)**
  - Higher Qb → thinner boundary layer → better kb → more clearance
- **Practical limit:** >450 mL/min risks collapsing the vascular access (fistula)

### Dialysate Flow Rate (Qd = 500 mL/min)

The rate at which fresh dialysate is pumped through the shell side.

- **Clinical range:** 300-800 mL/min
- **Typical prescription:** 500 mL/min
- **What it affects:**
  - Shell-side velocity → **dialysate-side mass transfer coefficient (kd)**
  - Higher Qd → better kd
- **Diminishing returns:** Since membrane resistance dominates (~70%), doubling Qd might only improve clearance by ~5-10%. The bottleneck is the membrane, not the dialysate film.

### Ultrafiltration Rate (Quf = 10 mL/min)

The net rate of water removal from the patient's blood.

- **Clinical range:** 5-15 mL/min
- **Purpose:** Remove excess body water accumulated between dialysis sessions
- **What it affects:**
  1. **Flow rates change along the fiber:** Qb(z) decreases, Qd(z) increases
  2. **Adds convective solute transport:** Solute gets dragged along with the water through the pores
  3. **Total fluid removed per session:** Quf × time = 10 × 240 min = 2,400 mL = 2.4 liters
- **Safety limit:** Too high Quf → blood pressure drops (hypotension), cramping, nausea

---

## 5. Calculated Parameters (Not Inputs — Derived from the Above)

These are not set directly but computed from the input parameters:

| Parameter | Formula | Typical Value | What It Represents |
|---|---|---|---|
| **D_water** (Urea) | kBT / (6πμr_s) | 1.18 × 10⁻⁹ m²/s | How fast urea diffuses in water |
| **D_m** (Urea) | D_water × (ε/τ) × (1-λ)² | ~2.5 × 10⁻¹⁰ m²/s | How fast urea diffuses through membrane |
| **kb** (Urea) | Sh × D / d_i | ~2.3 × 10⁻⁵ m/s | Blood-side film coefficient |
| **kd** (Urea) | Sh_d × D / d_h | ~7.8 × 10⁻⁵ m/s | Dialysate-side film coefficient |
| **Ko** (Urea) | 1/(1/kb + δ/Dm + 1/kd) | ~5.5 × 10⁻⁶ m/s | Overall mass transfer coefficient |
| **Sieving (S)** | (1-λ)²(2-(1-λ)²) | Urea: 0.97, Alb: 0 | Fraction convected through pores |
| **Total Area (A)** | n × π × d_i × L | 1.57 m² | Total membrane surface area |
| **μ_blood** | μ_water × (1 + 2.5×Hct) | 1.38 × 10⁻³ Pa·s | Blood viscosity |

---

## Quick Reference: What Happens When You Change Each Parameter

| Parameter | Increase → | Physical Reason |
|---|---|---|
| **Qb ↑** | K increases (more for small solutes) | Thinner blood boundary layer → higher kb |
| **Qd ↑** | K increases slightly | Thinner dialysate boundary layer → higher kd (but membrane dominates) |
| **Quf ↑** | K increases slightly + more water removed | Adds convective transport + extra water removal |
| **r_pore ↑** | Higher clearance for middle molecules, risk of albumin leak | More molecules fit through, steric hindrance decreases |
| **n_fibers ↑** | K increases | More membrane area = more transfer surface |
| **L ↑** | K increases | Longer contact time between blood and dialysate |
| **ε ↑** | K increases | More open channels for diffusion through membrane |
| **τ ↑** | K decreases | Longer, more tortuous path slows diffusion |
| **δ ↑** | K decreases | Thicker membrane = more resistance |
| **Hct ↑** | K decreases | Higher blood viscosity → thicker boundary layer |
