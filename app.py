"""
Separator Sizing Calculator
Based on API 12J (Specification for Oil and Gas Separators)
and GPSA Engineering Data Book
SI Units Throughout
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import math

# ─────────────────────────────────────────────────────────────────────────────
# PAGE CONFIG
# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="OpenSep",
    page_icon="⚙️",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─────────────────────────────────────────────────────────────────────────────
# CUSTOM CSS
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'IBM Plex Sans', sans-serif;
}

/* Dark industrial theme */
.stApp {
    background-color: #0d1117;
    color: #e6edf3;
}

/* Sidebar */
[data-testid="stSidebar"] {
    background-color: #161b22;
    border-right: 1px solid #30363d;
}

[data-testid="stSidebar"] .stMarkdown h2,
[data-testid="stSidebar"] .stMarkdown h3 {
    color: #f0a500;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.75rem;
    text-transform: uppercase;
    letter-spacing: 0.15em;
    border-bottom: 1px solid #30363d;
    padding-bottom: 0.3rem;
    margin-top: 1.2rem;
}

/* Cards */
.result-card {
    background: #161b22;
    border: 1px solid #30363d;
    border-left: 3px solid #f0a500;
    border-radius: 6px;
    padding: 1rem 1.25rem;
    margin: 0.4rem 0;
}
.result-card .label {
    font-size: 0.72rem;
    color: #8b949e;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    font-family: 'IBM Plex Mono', monospace;
}
.result-card .value {
    font-size: 1.5rem;
    font-weight: 700;
    color: #f0a500;
    font-family: 'IBM Plex Mono', monospace;
}
.result-card .unit {
    font-size: 0.85rem;
    color: #8b949e;
    margin-left: 0.3rem;
}

.warn-card {
    background: #1c1810;
    border: 1px solid #f0a500;
    border-radius: 6px;
    padding: 0.75rem 1rem;
    margin: 0.4rem 0;
    color: #f0a500;
    font-size: 0.85rem;
}
.ok-card {
    background: #0d1f12;
    border: 1px solid #238636;
    border-radius: 6px;
    padding: 0.75rem 1rem;
    margin: 0.4rem 0;
    color: #3fb950;
    font-size: 0.85rem;
}

/* Tabs */
.stTabs [data-baseweb="tab-list"] {
    background-color: #161b22;
    border-bottom: 1px solid #30363d;
    gap: 0;
}
.stTabs [data-baseweb="tab"] {
    background-color: transparent;
    color: #8b949e;
    border-bottom: 2px solid transparent;
    padding: 0.6rem 1.2rem;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.8rem;
    text-transform: uppercase;
    letter-spacing: 0.08em;
}
.stTabs [aria-selected="true"] {
    color: #f0a500 !important;
    border-bottom: 2px solid #f0a500 !important;
    background-color: transparent !important;
}

/* Headers */
h1 { 
    font-family: 'IBM Plex Mono', monospace !important; 
    color: #f0a500 !important;
    font-size: 1.6rem !important;
    letter-spacing: -0.02em;
}
h2, h3 { 
    font-family: 'IBM Plex Mono', monospace !important;
    color: #e6edf3 !important;
}

/* Inputs */
.stNumberInput input, .stSelectbox select, .stSlider {
    background-color: #21262d !important;
    border-color: #30363d !important;
    color: #e6edf3 !important;
}

/* Metric */
[data-testid="metric-container"] {
    background: #161b22;
    border: 1px solid #30363d;
    border-radius: 6px;
    padding: 0.8rem 1rem;
}
[data-testid="metric-container"] label {
    color: #8b949e !important;
    font-size: 0.72rem !important;
    font-family: 'IBM Plex Mono', monospace !important;
    text-transform: uppercase;
    letter-spacing: 0.1em;
}
[data-testid="metric-container"] [data-testid="metric-value"] {
    color: #f0a500 !important;
    font-family: 'IBM Plex Mono', monospace !important;
}

/* Divider */
hr { border-color: #30363d !important; }

/* Info box */
.stInfo { background-color: #1c2128; border-color: #388bfd; }

/* Section header pill */
.section-pill {
    display: inline-block;
    background: #f0a50020;
    border: 1px solid #f0a50060;
    color: #f0a500;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.68rem;
    text-transform: uppercase;
    letter-spacing: 0.15em;
    padding: 0.15rem 0.6rem;
    border-radius: 100px;
    margin-bottom: 0.6rem;
}
</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# CONSTANTS & HELPER FUNCTIONS (API 12J / GPSA)
# ─────────────────────────────────────────────────────────────────────────────

def sutherland_viscosity(T_K, species="gas"):
    """Sutherland's formula – not used directly but kept for reference."""
    pass

def gas_density(P_kPa, T_K, MW, Z=1.0):
    """ρ_g [kg/m³] from real gas law. R = 8314 J/(kmol·K)."""
    R = 8314.0
    return (P_kPa * 1000 * MW) / (Z * R * T_K)

def liquid_density_correction(rho_l_15, T_C):
    """Simple thermal expansion correction for crude (API 11.1 light version)."""
    # Δρ ≈ -0.65 kg/m³ per °C above 15 °C (rough average)
    return rho_l_15 - 0.65 * (T_C - 15)

def cunningham_correction(dp_m, P_kPa, T_K, MW_gas):
    """Cunningham slip correction (important for dp < 3 µm)."""
    # Mean free path λ [m]
    mu = gas_viscosity_estimate(T_K, MW_gas)
    rho_g = gas_density(P_kPa, T_K, MW_gas)
    lam = (mu / (0.499 * rho_g)) * math.sqrt(math.pi * MW_gas / (2 * 8314 * T_K))
    Kn = 2 * lam / dp_m
    Cc = 1 + Kn * (1.257 + 0.4 * math.exp(-1.1 / Kn))
    return Cc

def gas_viscosity_estimate(T_K, MW):
    """Estimate gas viscosity [Pa·s] using Lee-Kesler correlation (simplified)."""
    T_R = T_K / (MW ** (1/3))  # rough reduced temp
    # Chapman-Enskog simplified
    mu = (2.6693e-6 * math.sqrt(MW * T_K)) / (3.5 ** 2 * 1.22)
    return max(mu, 1e-6)

def terminal_velocity_stokes(dp_m, rho_l, rho_g, mu_g):
    """Stokes' law settling velocity [m/s] for liquid droplets in gas."""
    return (dp_m ** 2 * (rho_l - rho_g) * 9.81) / (18 * mu_g)

def terminal_velocity_cd(dp_m, rho_disp, rho_cont, mu_cont):
    """
    Iterative drag coefficient settling/rising velocity [m/s].
    Works for both gas/liquid AND liquid/liquid systems.
    Uses Schiller-Naumann CD correlation.
    API 12J Section 5 / GPSA Chapter 7.

    Args:
        dp_m    : droplet diameter [m]
        rho_disp: dispersed phase density [kg/m³]  (liquid drop or oil drop)
        rho_cont: continuous phase density [kg/m³] (gas, water, or oil)
        mu_cont : continuous phase viscosity [Pa·s]
    Returns:
        terminal velocity [m/s] (always positive; direction implied by caller)
    """
    g = 9.81
    delta_rho = abs(rho_disp - rho_cont)
    if delta_rho < 1e-3 or mu_cont <= 0 or dp_m <= 0:
        return 1e-6  # near-zero but non-zero to avoid 999 sentinel

    # Start from Stokes
    vt = (dp_m**2 * delta_rho * g) / (18.0 * mu_cont)
    for _ in range(200):
        Re = rho_cont * vt * dp_m / mu_cont
        Re = max(Re, 1e-10)
        if Re < 0.4:
            CD = 24.0 / Re
        elif Re < 500.0:
            CD = 24.0 / Re + 6.0 / (1.0 + math.sqrt(Re)) + 0.4
        else:
            CD = 0.44
        vt_new = math.sqrt((4.0 * delta_rho * g * dp_m) / (3.0 * CD * rho_cont))
        if abs(vt_new - vt) / (vt + 1e-12) < 1e-6:
            return vt_new
        vt = vt_new
    return vt

def souders_brown_K(T_K, P_kPa, service="gas_liquid"):
    """
    Souders-Brown coefficient K [m/s] per API 12J Table 1 / GPSA Fig 7-3.
    Corrects for pressure.
    """
    # Base K from GPSA (approximate from chart digitisation)
    P_psia = P_kPa * 0.145038
    if P_psia <= 100:
        K_base = 0.122
    elif P_psia <= 200:
        K_base = 0.116
    elif P_psia <= 400:
        K_base = 0.107
    elif P_psia <= 600:
        K_base = 0.098
    elif P_psia <= 900:
        K_base = 0.085
    elif P_psia <= 1200:
        K_base = 0.073
    else:
        K_base = 0.061
    # Demister pad factor
    return K_base  # m/s

def max_gas_velocity_horizontal(K, rho_g, rho_l):
    """Maximum allowable gas velocity [m/s] – Souders-Brown."""
    return K * math.sqrt((rho_l - rho_g) / rho_g)


# ─────────────────────────────────────────────────────────────────────────────
# INLET DEVICE PARAMETERS  (API 12J §6.3 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

INLET_DEVICES = {
    "Half-pipe / Diverter plate": {
        "rho_v2_limit": 1490,    # kg/(m·s²)  API 12J Eq. 6-1
        "K_derating":   1.00,    # no derating — baseline device
        "L_correction": 1.00,    # no length penalty
        "pre_sep_credit": 0.0,   # fraction of liquid pre-separated
        "desc": "Standard baseline device. API 12J ρv² ≤ 1490 kg/(m·s²).",
    },
    "Vane-type inlet (e.g. Schoepentoeter)": {
        "rho_v2_limit": 6000,
        "K_derating":   1.05,    # slight K improvement due to better flow distribution
        "L_correction": 0.85,    # 15% shorter vessel — better gas distribution
        "pre_sep_credit": 0.30,  # 30% of inlet liquid pre-separated
        "desc": "High-capacity vane inlet. Higher ρv² limit, 15% L reduction, 30% pre-sep credit.",
    },
    "Cyclonic inlet device": {
        "rho_v2_limit": 12000,
        "K_derating":   1.10,    # better K due to pre-cleaned gas
        "L_correction": 0.75,    # 25% shorter vessel
        "pre_sep_credit": 0.60,  # 60% pre-separation
        "desc": "Cyclonic inlet. Highest capacity, 25% L reduction, 60% liquid pre-separation.",
    },
    "Inlet nozzle only (no device)": {
        "rho_v2_limit": 1490,
        "K_derating":   0.80,    # 20% K derating per GPSA recommendation
        "L_correction": 1.30,    # 30% length penalty — poor flow distribution
        "pre_sep_credit": 0.0,
        "desc": "No inlet device. K derated 20%, vessel length increased 30% per GPSA.",
    },
    "Perforated pipe / Spreader": {
        "rho_v2_limit": 2500,
        "K_derating":   0.95,
        "L_correction": 0.95,
        "pre_sep_credit": 0.10,
        "desc": "Perforated pipe distributor. Modest improvements over bare nozzle.",
    },
}

def get_inlet_device_params(device_name):
    """Return inlet device parameter dict."""
    return INLET_DEVICES.get(device_name, INLET_DEVICES["Half-pipe / Diverter plate"])


# ─────────────────────────────────────────────────────────────────────────────
# MIST EXTRACTOR DATA & SIZING  (API 12J §8 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

MIST_EXTRACTORS = {
    "None (gravity separation only)": {
        "K_factor":        None,   # uses Souders-Brown from pressure table, derated
        "K_derating":      0.75,   # 25% K derating without mist extractor
        "dp_micron":       None,   # no defined cut size
        "thickness_m":     0.0,    # no height contribution
        "dP_clean_Pa":     0.0,
        "dP_loaded_Pa":    0.0,
        "wire_diameter_mm":None,
        "void_fraction":   None,
        "specific_area_m2m3": None,
        "desc": "No mist extractor. K derated 25%. Suitable only where outlet gas quality is not critical.",
    },
    "Wire mesh pad (standard)": {
        "K_factor":        0.107,  # m/s — GPSA Fig. 7-3 at moderate pressure
        "K_derating":      1.00,
        "dp_micron":       10.0,   # typical cut size [µm]
        "thickness_m":     0.15,   # standard pad depth
        "dP_clean_Pa":     125.0,  # Pa clean pressure drop
        "dP_loaded_Pa":    375.0,  # Pa loaded (design)
        "wire_diameter_mm":0.28,   # typical wire diameter
        "void_fraction":   0.97,   # typical
        "specific_area_m2m3": 200, # m²/m³
        "desc": "Standard knitted wire mesh. dp_cut ≈ 10 µm. GPSA K=0.107 m/s. API 12J §8.3.",
    },
    "Wire mesh pad (high capacity)": {
        "K_factor":        0.122,
        "K_derating":      1.00,
        "dp_micron":       20.0,
        "thickness_m":     0.20,
        "dP_clean_Pa":     100.0,
        "dP_loaded_Pa":    300.0,
        "wire_diameter_mm":0.38,
        "void_fraction":   0.98,
        "specific_area_m2m3": 130,
        "desc": "High-capacity coarser mesh. Higher K, lower dP, larger cut size. For high gas rates.",
    },
    "Vane pack (chevron)": {
        "K_factor":        0.122,
        "K_derating":      1.00,
        "dp_micron":       30.0,
        "thickness_m":     0.30,
        "dP_clean_Pa":     200.0,
        "dP_loaded_Pa":    500.0,
        "wire_diameter_mm":None,
        "void_fraction":   0.90,
        "specific_area_m2m3": 250,
        "desc": "Vane/chevron pack. Higher dP, handles slugging better than mesh. dp_cut ≈ 30 µm.",
    },
    "Axial cyclone bundle": {
        "K_factor":        0.150,
        "K_derating":      1.00,
        "dp_micron":       5.0,
        "thickness_m":     0.50,
        "dP_clean_Pa":     500.0,
        "dP_loaded_Pa":    1200.0,
        "wire_diameter_mm":None,
        "void_fraction":   None,
        "specific_area_m2m3": None,
        "desc": "Axial cyclone tubes. Best efficiency (dp_cut ≈ 5 µm), highest K and dP. Offshore / high-spec.",
    },
    "Agglomeration filter": {
        "K_factor":        0.046,
        "K_derating":      1.00,
        "dp_micron":       3.0,
        "thickness_m":     0.25,
        "dP_clean_Pa":     250.0,
        "dP_loaded_Pa":    700.0,
        "wire_diameter_mm":None,
        "void_fraction":   0.92,
        "specific_area_m2m3": 500,
        "desc": "Agglomeration / coalescing filter. Sub-micron capability. Lowest K; use for critical duty.",
    },
}

def size_mist_extractor(me_type, Qg_m3s, rho_g, rho_l, D_vessel,
                         orientation="Horizontal", K_sb=None):
    """
    Size the mist extractor element per API 12J §8 / GPSA Ch. 7.

    Returns a dict with:
      - required area [m²]
      - actual area [m²] (based on vessel cross-section or annular area)
      - gas velocity through extractor [m/s]
      - design K [m/s]
      - allowable velocity [m/s]
      - utilisation [%]
      - pressure drop estimates [Pa]
      - efficiency note
      - all ME properties
    """
    me = MIST_EXTRACTORS.get(me_type, MIST_EXTRACTORS["Wire mesh pad (standard)"])

    if me_type == "None (gravity separation only)":
        return {
            "me_type": me_type,
            "desc": me["desc"],
            "K_design": None,
            "v_allow_ms": None,
            "A_required_m2": None,
            "A_actual_m2": None,
            "v_me_ms": None,
            "utilisation_pct": None,
            "dP_clean_Pa": 0.0,
            "dP_loaded_Pa": 0.0,
            "thickness_m": 0.0,
            "dp_cut_micron": None,
            "wire_dia_mm": None,
            "void_fraction": None,
            "specific_area": None,
            "K_derating": me["K_derating"],
            "no_me": True,
        }

    # Design K — use ME table value; if K_sb provided use whichever is lower
    K_design = me["K_factor"]
    if K_sb is not None:
        K_design = min(K_design, K_sb)

    # Allowable velocity through extractor
    v_allow = K_design * math.sqrt((rho_l - rho_g) / rho_g)

    # Required cross-sectional area
    A_required = Qg_m3s / v_allow

    # Actual available area
    A_vessel = math.pi * D_vessel**2 / 4
    if orientation == "Vertical":
        # Full vessel cross-section (minus ~10% for support ring)
        A_actual = 0.90 * A_vessel
    else:
        # Horizontal: ME sits in gas space (top 40-50%)
        A_actual = 0.45 * A_vessel

    # Actual velocity through ME
    v_me = Qg_m3s / A_actual

    # Utilisation
    utilisation = (v_me / v_allow) * 100.0

    # If actual area < required — ME undersized for this vessel diameter
    undersized = A_actual < A_required

    return {
        "me_type": me_type,
        "desc": me["desc"],
        "K_design": round(K_design, 4),
        "v_allow_ms": round(v_allow, 4),
        "A_required_m2": round(A_required, 4),
        "A_actual_m2": round(A_actual, 4),
        "v_me_ms": round(v_me, 4),
        "utilisation_pct": round(utilisation, 1),
        "dP_clean_Pa": me["dP_clean_Pa"],
        "dP_loaded_Pa": me["dP_loaded_Pa"],
        "thickness_m": me["thickness_m"],
        "dp_cut_micron": me["dp_micron"],
        "wire_dia_mm": me["wire_diameter_mm"],
        "void_fraction": me["void_fraction"],
        "specific_area": me["specific_area_m2m3"],
        "K_derating": me["K_derating"],
        "undersized": undersized,
        "no_me": False,
    }

def retention_time_liquid(API_gravity, GOR, phase="two"):
    """
    Recommended liquid retention time [s] per API 12J Table 2.
    API_gravity: oil API gravity
    GOR: gas-oil ratio [m³/m³] (std)
    """
    # Convert GOR to scf/bbl
    GOR_scfbbl = GOR * 5.6146  # 1 m³/m³ ≈ 5.6146 scf/bbl? No, GOR m³(gas)/m³(liq)
    # API 12J Table 2: retention times
    if API_gravity > 35:
        t_r = 3 * 60  # 3 min
    elif API_gravity > 25:
        t_r = 5 * 60  # 5 min
    else:
        t_r = 10 * 60  # 10 min
    if phase == "three":
        t_r = max(t_r, 10 * 60)
    return t_r  # seconds

def water_retention_time(T_C):
    """Water retention time [s] per API 12J."""
    if T_C > 60:
        return 5 * 60
    elif T_C > 40:
        return 10 * 60
    else:
        return 20 * 60

# ─────────────────────────────────────────────────────────────────────────────
# HORIZONTAL SEPARATOR SIZING  (API 12J Sec. 5.3 + GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

def size_horizontal_two_phase(
    Qg_m3s,   # gas volumetric flow [m³/s] at operating conditions
    Ql_m3s,   # liquid flow [m³/s] at operating conditions
    rho_g,    # gas density [kg/m³]
    rho_l,    # liquid density [kg/m³]
    mu_g,     # gas viscosity [Pa·s]
    dp_mic,   # design droplet diameter [µm]
    t_ret_s,  # liquid retention time [s]
    T_K,      # temperature [K]
    P_kPa,    # pressure [kPa]
    K_override=None,  # optional manual K [m/s]
    inlet_device="Half-pipe / Diverter plate",
    me_K_derating=1.0,
):
    """
    Returns dict of sizing results for a horizontal two-phase separator.
    Methodology: API 12J Section 5.3 / GPSA Engineering Data Book Chapter 7.
    """
    dp_m = dp_mic * 1e-6
    g = 9.81
    dev = get_inlet_device_params(inlet_device)

    # 1. Terminal settling velocity of liquid droplets in gas
    vt = terminal_velocity_cd(dp_m, rho_l, rho_g, mu_g)

    # 2. Souders-Brown K — apply inlet device + ME derating
    K_base = K_override if K_override is not None else souders_brown_K(T_K, P_kPa)
    K = K_base * dev["K_derating"] * me_K_derating
    v_max_SB = max_gas_velocity_horizontal(K, rho_g, rho_l)

    # 3. Liquid volume requirement
    V_liq = Ql_m3s * t_ret_s   # [m³]

    # 4. Try candidate diameters (0.3 m to 5.0 m, step 0.05 m)
    results = []
    for D in np.arange(0.3, 5.05, 0.05):
        # Assume 50% liquid level for two-phase (API 12J recommendation)
        # Cross-section areas
        A_total = math.pi * D ** 2 / 4
        # Gas occupies top half, liquid bottom half (initial assumption)
        frac_liq = 0.50
        A_liq = frac_liq * A_total
        A_gas = (1 - frac_liq) * A_total

        # Effective gas velocity
        v_gas = Qg_m3s / A_gas

        if v_gas > v_max_SB:  # must not exceed Souders-Brown limit
            continue

        # Droplet settling check: L_eff ≥ v_gas * D / vt
        # From API 12J: L_eff / D ≥ v_gas / vt  (for horizontal flow)
        L_min_gas = (v_gas * D) / vt

        # Liquid retention length
        # V_liq = A_liq * L_eff  → L_eff = V_liq / A_liq
        L_min_liq = V_liq / A_liq

        # Without ME, gas settling length is penalised: droplets must
        # travel further to gravity-settle to the liquid surface
        L_gas_me_factor = 1.0 if me_K_derating >= 1.0 else (1.0 / me_K_derating)
        L_eff = max(L_min_gas * L_gas_me_factor, L_min_liq) * dev["L_correction"]

        # L/D ratio (API 12J: 3 ≤ L/D ≤ 5 preferred, allow up to 6)
        LD = L_eff / D

        if LD < 1.5:
            continue  # unrealistically short

        # Seam-to-seam length = L_eff + 2*(0.15*D) heads (GPSA rule)
        L_ss = L_eff + 0.3 * D

        # Check droplet removal efficiency (simplified)
        efficiency = min(1.0, vt * L_eff / (v_gas * D))

        results.append({
            "D_m": round(D, 2),
            "L_eff_m": round(L_eff, 3),
            "L_ss_m": round(L_ss, 3),
            "LD": round(LD, 2),
            "v_gas_ms": round(v_gas, 3),
            "v_max_ms": round(v_max_SB, 3),
            "vt_ms": round(vt, 5),
            "efficiency": round(efficiency * 100, 1),
            "V_liq_m3": round(V_liq, 3),
            "K_factor": round(K, 4),
            "inlet_device": inlet_device,
            "K_derating": round(dev["K_derating"], 3),
            "L_correction": round(dev["L_correction"], 3),
        })

    if not results:
        return None, None

    # Select preferred: smallest vessel satisfying 3 ≤ L/D ≤ 5
    preferred = [r for r in results if 3.0 <= r["LD"] <= 5.0]
    if not preferred:
        preferred = [r for r in results if r["LD"] <= 6.0]
    if not preferred:
        preferred = results

    best = preferred[0]
    return best, pd.DataFrame(results)


def size_horizontal_three_phase(
    Qg_m3s, Qo_m3s, Qw_m3s,
    rho_g, rho_o, rho_w,
    mu_g, mu_o, mu_w,
    dp_o_mic, dp_w_mic,
    t_ret_o_s, t_ret_w_s,
    T_K, P_kPa,
    K_override=None,
    inlet_device="Half-pipe / Diverter plate",
    me_K_derating=1.0,
):
    """
    Horizontal three-phase separator sizing.
    API 12J Section 5.4 / GPSA Chapter 7.
    """
    dp_o = dp_o_mic * 1e-6
    dp_w = dp_w_mic * 1e-6
    g = 9.81

    vt_gas = terminal_velocity_cd(dp_o, rho_o, rho_g, mu_g)   # oil drop in gas
    vt_o_in_w = terminal_velocity_cd(dp_o, rho_o, rho_w, mu_w) # oil drop rise in water
    vt_w_in_o = terminal_velocity_cd(dp_w, rho_w, rho_o, mu_o) # water drop fall in oil

    dev = get_inlet_device_params(inlet_device)
    K_base = K_override if K_override is not None else souders_brown_K(T_K, P_kPa)
    K = K_base * dev["K_derating"] * me_K_derating
    v_max_SB = max_gas_velocity_horizontal(K, rho_g, rho_o)

    # Pre-separation credit: reduce liquid load entering vessel
    pre_sep = dev["pre_sep_credit"]
    Qo_m3s_eff = Qo_m3s * (1.0 - pre_sep)
    Qw_m3s_eff = Qw_m3s * (1.0 - pre_sep)

    # Volumes
    V_o = Qo_m3s_eff * t_ret_o_s
    V_w = Qw_m3s_eff * t_ret_w_s
    V_liq = V_o + V_w

    results = []
    for D in np.arange(0.3, 5.55, 0.05):
        A_total = math.pi * D ** 2 / 4

        # Phase fraction assignments (API 12J Fig 5-2 guide)
        # Gas: top 40%, Oil: middle 35%, Water: bottom 25%
        f_gas = 0.40
        f_oil = 0.35
        f_wat = 0.25

        A_gas = f_gas * A_total
        A_oil = f_oil * A_total
        A_wat = f_wat * A_total

        v_gas = Qg_m3s / A_gas
        if v_gas > v_max_SB:
            continue

        # Gas section: oil droplet settling
        L_gas = (v_gas * D) / vt_gas

        # Oil section: water droplet settling (API 12J Eq. for coalescing)
        v_oil = Qo_m3s_eff / A_oil
        L_oil_w = (v_oil * D) / vt_w_in_o

        # Water section: oil droplet rise
        v_wat = Qw_m3s_eff / A_wat
        L_wat_o = (v_wat * D) / vt_o_in_w

        # Retention lengths
        L_ret_o = V_o / A_oil
        L_ret_w = V_w / A_wat

        L_gas_me_factor = 1.0 if me_K_derating >= 1.0 else (1.0 / me_K_derating)
        L_eff = max(L_gas * L_gas_me_factor, L_oil_w, L_wat_o, L_ret_o, L_ret_w) * dev["L_correction"]
        LD = L_eff / D

        if LD < 1.5:
            continue

        L_ss = L_eff + 0.3 * D

        results.append({
            "D_m": round(D, 2),
            "L_eff_m": round(L_eff, 3),
            "L_ss_m": round(L_ss, 3),
            "LD": round(LD, 2),
            "v_gas_ms": round(v_gas, 4),
            "v_max_ms": round(v_max_SB, 3),
            "L_gas_ctrl_m": round(L_gas, 3),
            "L_oil_ctrl_m": round(L_oil_w, 3),
            "L_wat_ctrl_m": round(L_wat_o, 3),
            "V_oil_m3": round(V_o, 3),
            "V_wat_m3": round(V_w, 3),
            "K_factor": round(K, 4),
            "vt_ms": round(vt_gas, 5),
            "inlet_device": inlet_device,
            "K_derating": round(dev["K_derating"], 3),
            "L_correction": round(dev["L_correction"], 3),
            "pre_sep_credit": round(pre_sep * 100, 0),
        })

    if not results:
        return None, None

    preferred = [r for r in results if 3.0 <= r["LD"] <= 5.0]
    if not preferred:
        preferred = [r for r in results if r["LD"] <= 6.0]
    if not preferred:
        preferred = results

    best = preferred[0]
    return best, pd.DataFrame(results)


# ─────────────────────────────────────────────────────────────────────────────
# VERTICAL SEPARATOR SIZING  (API 12J Sec. 5.2 + GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

def size_vertical_two_phase(
    Qg_m3s, Ql_m3s,
    rho_g, rho_l,
    mu_g, dp_mic,
    t_ret_s, T_K, P_kPa,
    K_override=None,
    inlet_device="Half-pipe / Diverter plate",
    me_thickness=0.15,
    me_K_derating=1.0,
):
    """
    Vertical two-phase separator sizing.
    API 12J Section 5.2 / GPSA Chapter 7.
    """
    dp_m = dp_mic * 1e-6
    dev = get_inlet_device_params(inlet_device)
    vt = terminal_velocity_cd(dp_m, rho_l, rho_g, mu_g)
    K_base = K_override if K_override is not None else souders_brown_K(T_K, P_kPa)
    K = K_base * dev["K_derating"] * me_K_derating
    v_max_SB = K * math.sqrt((rho_l - rho_g) / rho_g)

    # For vertical: design gas velocity ≤ 75% * vt (GPSA Eq. 7-4) or Souders-Brown
    v_design = min(v_max_SB, 0.75 * vt)

    # Required cross-sectional area for gas
    A_req = Qg_m3s / v_design
    D_min_gas = math.sqrt(4 * A_req / math.pi)

    # Liquid sump volume
    V_liq = Ql_m3s * t_ret_s

    results = []
    for D in np.arange(0.3, 5.05, 0.05):
        A = math.pi * D ** 2 / 4
        v_gas = Qg_m3s / A

        if v_gas > v_max_SB:
            continue

        # Liquid sump height: H_liq = V_liq / A  (cylindrical sump)
        H_liq = V_liq / A

        # Gas disengagement height: API 12J §5.2
        # Without ME: ≥ 0.9 m (gravity settling must do more work)
        # With ME   : ≥ 0.6 m
        H_gas_min = 0.6 if me_thickness > 0 else 0.9
        H_gas = max(H_gas_min, D * 0.5)

        # Mist extractor height (from ME selection)
        H_mist = me_thickness

        # Inlet nozzle zone: 0.3 m
        H_inlet = 0.3

        H_total = H_liq + H_gas + H_mist + H_inlet

        # Slenderness ratio H/D
        HD = H_total / D

        results.append({
            "D_m": round(D, 2),
            "H_liq_m": round(H_liq, 3),
            "H_gas_m": round(H_gas, 3),
            "H_mist_m": round(H_mist, 3),
            "H_total_m": round(H_total, 3),
            "HD": round(HD, 2),
            "v_gas_ms": round(v_gas, 4),
            "v_max_ms": round(v_max_SB, 4),
            "vt_ms": round(vt, 5),
            "V_liq_m3": round(V_liq, 3),
            "K_factor": round(K, 4),
            "inlet_device": inlet_device,
            "K_derating": round(dev["K_derating"], 3),
        })

    if not results:
        return None, None

    # Prefer H/D between 1.5 and 4 (API 12J / GPSA recommendation)
    preferred = [r for r in results if 1.5 <= r["HD"] <= 4.0]
    if not preferred:
        preferred = results

    best = preferred[0]
    return best, pd.DataFrame(results)


def size_vertical_three_phase(
    Qg_m3s, Qo_m3s, Qw_m3s,
    rho_g, rho_o, rho_w,
    mu_g, mu_o, mu_w,
    dp_o_mic, dp_w_mic,
    t_ret_o_s, t_ret_w_s,
    T_K, P_kPa,
    K_override=None,
    inlet_device="Half-pipe / Diverter plate",
    me_thickness=0.15,
    me_K_derating=1.0,
):
    """
    Vertical three-phase separator sizing.
    API 12J Section 5.4 / GPSA Chapter 7.
    """
    dp_o = dp_o_mic * 1e-6
    dp_w = dp_w_mic * 1e-6
    dev = get_inlet_device_params(inlet_device)

    vt_gas = terminal_velocity_cd(dp_o, rho_o, rho_g, mu_g)
    vt_o_w = terminal_velocity_cd(dp_o, rho_o, rho_w, mu_w)  # oil rise in water

    K_base = K_override if K_override is not None else souders_brown_K(T_K, P_kPa)
    K = K_base * dev["K_derating"] * me_K_derating
    v_max_SB = K * math.sqrt((rho_o - rho_g) / rho_g)
    v_design = min(v_max_SB, 0.75 * vt_gas)

    pre_sep = dev["pre_sep_credit"]
    V_o = Qo_m3s * (1.0 - pre_sep) * t_ret_o_s
    V_w = Qw_m3s * (1.0 - pre_sep) * t_ret_w_s

    results = []
    for D in np.arange(0.3, 5.55, 0.05):
        A = math.pi * D ** 2 / 4
        v_gas = Qg_m3s / A
        if v_gas > v_max_SB:
            continue

        # Oil pad thickness (coalescence zone) – API 12J
        # H_oil = V_o / A
        H_oil = V_o / A

        # Water leg height
        H_water = V_w / A

        # Gas disengagement zone
        H_gas_min = 0.6 if me_thickness > 0 else 0.9
        H_gas = max(H_gas_min, D * 0.5)
        H_mist = me_thickness
        H_inlet = 0.3

        H_total = H_water + H_oil + H_gas + H_mist + H_inlet
        HD = H_total / D

        results.append({
            "D_m": round(D, 2),
            "H_oil_m": round(H_oil, 3),
            "H_water_m": round(H_water, 3),
            "H_gas_m": round(H_gas, 3),
            "H_mist_m": round(H_mist, 3),
            "H_total_m": round(H_total, 3),
            "HD": round(HD, 2),
            "v_gas_ms": round(v_gas, 4),
            "v_max_ms": round(v_max_SB, 4),
            "vt_ms": round(vt_gas, 5),
            "vt_gas_ms": round(vt_gas, 5),
            "vt_o_w_ms": round(vt_o_w, 5),
            "V_oil_m3": round(V_o, 3),
            "V_wat_m3": round(V_w, 3),
            "K_factor": round(K, 4),
        })

    if not results:
        return None, None

    preferred = [r for r in results if 1.5 <= r["HD"] <= 4.0]
    if not preferred:
        preferred = results

    best = preferred[0]
    return best, pd.DataFrame(results)



# ─────────────────────────────────────────────────────────────────────────────
# NOZZLE SIZING & POSITIONS  (API 12J Sec. 6 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

# Standard nozzle OD sizes [m] per ASME B36.10 / API 6D (DN 50 → DN 600)
STANDARD_NOZZLE_DN = [
    0.0508, 0.0635, 0.0762, 0.1016, 0.1524, 0.2032,
    0.2540, 0.3048, 0.3556, 0.4064, 0.5080, 0.6096
]

def next_standard_nozzle(d_calc):
    """Return next standard DN nozzle ID [m] >= d_calc."""
    for d in STANDARD_NOZZLE_DN:
        if d >= d_calc:
            return d
    return STANDARD_NOZZLE_DN[-1]

def nozzle_diameter(Q_m3s, rho, nozzle_type="liquid", rho_v2_limit=1490):
    """
    Size a nozzle per API 12J Section 6 momentum / velocity criteria.

    Inlet  : ρv² ≤ 1490 kg/(m·s²)  [API 12J §6.3, Eq. 6-1]
    Gas out: v ≤ 18 m/s             [API 12J §6.4 guideline]
    Liq out: v ≤  3 m/s             [API 12J §6.5 guideline]
    Relief  : v ≤ 18 m/s
    """
    if Q_m3s <= 0 or rho <= 0:
        return {"d_calc_m": 0.0, "d_sel_m": STANDARD_NOZZLE_DN[0],
                "velocity_ms": 0.0, "momentum": 0.0, "criterion": "—"}

    if nozzle_type == "inlet":
        # API 12J Eq. 6-1:  ρv² ≤ rho_v2_limit  →  A ≥ Q * sqrt(ρ/limit)
        # d = sqrt(4A/π)
        rho_v2_limit = float(rho_v2_limit)  # from inlet device selection
        A_min = Q_m3s * math.sqrt(rho / rho_v2_limit)
        d_calc = math.sqrt(4 * A_min / math.pi)
        d_sel = next_standard_nozzle(d_calc)
        A_sel = math.pi * d_sel**2 / 4
        v = Q_m3s / A_sel
        rho_v2 = rho * v**2
        return {"d_calc_m": d_calc, "d_sel_m": d_sel,
                "velocity_ms": v, "momentum": rho_v2,
                "criterion": f"ρv²={rho_v2:.0f} ≤ 1490 kg/(m·s²)"}

    elif nozzle_type == "gas_outlet":
        v_limit = 18.0
        A_min = Q_m3s / v_limit
        d_calc = math.sqrt(4 * A_min / math.pi)
        d_sel = next_standard_nozzle(d_calc)
        A_sel = math.pi * d_sel**2 / 4
        v = Q_m3s / A_sel
        return {"d_calc_m": d_calc, "d_sel_m": d_sel,
                "velocity_ms": v, "momentum": rho * v**2,
                "criterion": f"v={v:.2f} ≤ 18 m/s"}

    elif nozzle_type == "liquid_outlet":
        v_limit = 3.0
        A_min = Q_m3s / v_limit
        d_calc = math.sqrt(4 * A_min / math.pi)
        d_sel = next_standard_nozzle(d_calc)
        A_sel = math.pi * d_sel**2 / 4
        v = Q_m3s / A_sel
        return {"d_calc_m": d_calc, "d_sel_m": d_sel,
                "velocity_ms": v, "momentum": rho * v**2,
                "criterion": f"v={v:.2f} ≤ 3 m/s"}

    elif nozzle_type == "water_outlet":
        v_limit = 3.0
        A_min = Q_m3s / v_limit
        d_calc = math.sqrt(4 * A_min / math.pi)
        d_sel = next_standard_nozzle(d_calc)
        A_sel = math.pi * d_sel**2 / 4
        v = Q_m3s / A_sel
        return {"d_calc_m": d_calc, "d_sel_m": d_sel,
                "velocity_ms": v, "momentum": rho * v**2,
                "criterion": f"v={v:.2f} ≤ 3 m/s"}

    elif nozzle_type == "relief":
        v_limit = 18.0
        A_min = Q_m3s / v_limit
        d_calc = math.sqrt(4 * A_min / math.pi)
        d_sel = next_standard_nozzle(d_calc)
        A_sel = math.pi * d_sel**2 / 4
        v = Q_m3s / A_sel
        return {"d_calc_m": d_calc, "d_sel_m": d_sel,
                "velocity_ms": v, "momentum": rho * v**2,
                "criterion": f"v={v:.2f} ≤ 18 m/s"}

    return {"d_calc_m": 0.0, "d_sel_m": STANDARD_NOZZLE_DN[0],
            "velocity_ms": 0.0, "momentum": 0.0, "criterion": "—"}


def nozzle_positions_horizontal(D, L_ss, best, is_three, n_inlet=1):
    """
    Compute nozzle centreline positions for a horizontal separator.
    Reference: API 12J §6 and industry standard layout practice.

    Coordinate system:
      x : axial (0 = inlet head seam, L_ss = outlet head seam)
      y : elevation above vessel bottom (0 = bottom, D = top)

    Returns dict of nozzle positions {name: (x, y, description)}.
    """
    L_eff = best.get("L_eff_m", L_ss - 0.15*D)
    head_h = 0.15 * D   # head tangent length (approx)

    # Liquid levels from sizing result
    if is_three:
        # Water level = 25% of D, oil/water interface, oil level = 60% of D
        y_water = D * 0.25
        y_oil   = D * 0.60
        y_hll   = D * 0.62   # high liquid level
        y_lll   = D * 0.10   # low liquid level
    else:
        y_liq   = D * 0.50   # normal liquid level (50% fill)
        y_hll   = D * 0.55
        y_lll   = D * 0.15

    # ── Positions ─────────────────────────────────────────────────────────────
    pos = {}

    # Inlet nozzle: on inlet head, centred at D/2 height (slightly above NLL)
    x_inlet = head_h
    y_inlet = D * 0.65   # above liquid level to allow gas flash
    pos["Inlet"] = (x_inlet, y_inlet, f"x={x_inlet:.3f} m from inlet seam, y={y_inlet:.3f} m from bottom")

    # Gas outlet: top of vessel, near outlet end
    x_gas_out = L_ss - head_h
    y_gas_out = D * 0.92
    pos["Gas Outlet"] = (x_gas_out, y_gas_out, f"x={x_gas_out:.3f} m, y={y_gas_out:.3f} m (top)")

    if is_three:
        # Oil outlet: outlet end, at oil level (weir controlled)
        x_oil_out = L_ss - head_h
        y_oil_out = y_oil
        pos["Oil Outlet"] = (x_oil_out, y_oil_out, f"x={x_oil_out:.3f} m, y={y_oil_out:.3f} m")

        # Water outlet: outlet end, near bottom
        x_wat_out = L_ss - head_h
        y_wat_out = D * 0.08
        pos["Water Outlet"] = (x_wat_out, y_wat_out, f"x={x_wat_out:.3f} m, y={y_wat_out:.3f} m")
    else:
        # Liquid outlet: outlet end, at liquid level
        x_liq_out = L_ss - head_h
        y_liq_out = D * 0.15
        pos["Liquid Outlet"] = (x_liq_out, y_liq_out, f"x={x_liq_out:.3f} m, y={y_liq_out:.3f} m")

    # Drain: bottom centre
    pos["Drain"] = (L_ss * 0.5, 0.0, f"x={L_ss*0.5:.3f} m, bottom of vessel")

    # Relief / vent: top near inlet
    x_relief = head_h + 0.3
    pos["Relief/Vent"] = (x_relief, D * 0.95, f"x={x_relief:.3f} m, top of vessel")

    # Level instrument (LG): inlet end
    pos["Level Gauge"] = (head_h + 0.2, D * 0.50, f"x={head_h+0.2:.3f} m, at NLL")

    return pos


def weir_position_horizontal(D, L_ss, best, is_three,
                              Qo_m3s=0, Qw_m3s=0,
                              rho_o=800, rho_w=1050):
    """
    Compute weir height and position for a horizontal three-phase separator.
    API 12J §5.4.3 / GPSA §7.

    Weir height h_w controls the oil/water interface level.
    h_w is set so the water retention volume is maintained.

    Returns dict with weir x-position, height, and oil pad thickness.
    """
    if not is_three:
        return None

    A_total = math.pi * D**2 / 4
    f_wat = 0.25
    A_wat = f_wat * A_total
    V_wat = best.get("V_wat_m3", 0)

    # Weir is located at the outlet end of the liquid section
    # x_weir ≈ L_eff (separates liquid compartment from outlet)
    L_eff = best.get("L_eff_m", L_ss * 0.85)
    x_weir = L_eff

    # Weir height = water level = 25% of D (as set in sizing)
    h_weir = D * 0.25   # m above vessel bottom

    # Oil pad above the weir
    h_oil_pad = D * 0.60 - h_weir   # oil/water interface to oil outlet level

    # Oil overflow weir: oil level = 60% D, so weir top for oil = 60% D
    h_oil_weir = D * 0.60

    return {
        "x_weir_m": round(x_weir, 3),
        "h_water_weir_m": round(h_weir, 3),
        "h_oil_weir_m": round(h_oil_weir, 3),
        "oil_pad_thickness_m": round(h_oil_pad, 3),
        "desc": (
            f"Water weir at x={x_weir:.3f} m, height={h_weir:.3f} m.  "
            f"Oil overflow weir height={h_oil_weir:.3f} m.  "
            f"Oil pad thickness={h_oil_pad:.3f} m."
        )
    }


def nozzle_positions_vertical(D, H_total, best, is_three):
    """
    Nozzle positions for a vertical separator.
    Elevation z measured from bottom tangent line.
    """
    H_liq  = best.get("H_liq_m",  best.get("H_oil_m", 0) + best.get("H_water_m", 0))
    H_gas  = best.get("H_gas_m",  max(0.6, D * 0.5))
    H_mist = 0.15
    H_inlet_zone = 0.3

    z_bottom    = 0.0
    z_nll       = H_liq * 0.5     # normal liquid level
    z_hll       = H_liq * 0.85
    z_lll       = H_liq * 0.10
    z_inlet     = H_liq + H_inlet_zone * 0.5   # mid of inlet zone
    z_mist_bot  = H_liq + H_gas
    z_gas_out   = H_total - 0.1

    pos = {}
    pos["Inlet"]       = (D*0.5, z_inlet,    f"z={z_inlet:.3f} m, side nozzle")
    pos["Gas Outlet"]  = (D*0.5, z_gas_out,  f"z={z_gas_out:.3f} m (top)")
    pos["Drain"]       = (D*0.5, 0.0,         "z=0.0 m, bottom")
    pos["Relief/Vent"] = (D*0.5, z_gas_out,   f"z={z_gas_out:.3f} m, top")
    pos["Level Gauge"] = (D*0.5, z_nll,        f"z={z_nll:.3f} m, at NLL")

    if is_three:
        z_oil_out   = best.get("H_water_m", 0) + best.get("H_oil_m", 0) * 0.5
        z_water_out = best.get("H_water_m", 0) * 0.15
        pos["Oil Outlet"]   = (D*0.5, z_oil_out,   f"z={z_oil_out:.3f} m")
        pos["Water Outlet"] = (D*0.5, z_water_out, f"z={z_water_out:.3f} m")
    else:
        z_liq_out = H_liq * 0.1
        pos["Liquid Outlet"] = (D*0.5, z_liq_out, f"z={z_liq_out:.3f} m")

    return pos

# ─────────────────────────────────────────────────────────────────────────────
# PLOT HELPERS
# ─────────────────────────────────────────────────────────────────────────────

PLOTLY_THEME = dict(
    paper_bgcolor="#0d1117",
    plot_bgcolor="#161b22",
    font=dict(family="IBM Plex Mono", color="#8b949e", size=11),
    title_font=dict(color="#e6edf3", size=13),
)

# Axis defaults for charts that use a standard grid
_AXIS_DEFAULTS = dict(
    xaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
    yaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
)

def plot_LD_curve(df, best, xlabel="D [m]", ylabel="L/D or H/D", title="Sizing Envelope"):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["D_m"], y=df["LD"] if "LD" in df else df["HD"],
        mode="lines", name="L/D ratio",
        line=dict(color="#388bfd", width=2)
    ))
    # Preferred band
    fig.add_hrect(y0=3, y1=5, fillcolor="rgba(240,165,0,0.08)",
                  line=dict(color="rgba(240,165,0,0.25)"), annotation_text="API 12J preferred",
                  annotation_font=dict(color="#f0a500", size=10))
    fig.add_vline(x=best["D_m"], line=dict(color="#f0a500", dash="dash", width=1.5),
                  annotation_text=f"Selected D={best['D_m']} m",
                  annotation_font=dict(color="#f0a500", size=10))
    fig.update_layout(
        **PLOTLY_THEME, **_AXIS_DEFAULTS,
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        height=320, margin=dict(l=50, r=20, t=40, b=40),
    )
    return fig

def plot_vessel_schematic_horizontal(D, L_ss, n_phases=2):
    """Simple 2-D cross-section schematic."""
    fig = go.Figure()
    # Shell outline
    fig.add_shape(type="rect", x0=0, y0=0, x1=L_ss, y1=D,
                  line=dict(color="#f0a500", width=2), fillcolor="#21262d")
    # Liquid level(s)
    if n_phases == 2:
        fig.add_shape(type="rect", x0=0.05*L_ss, y0=0, x1=0.95*L_ss, y1=D*0.50,
                      line=dict(width=0), fillcolor="rgba(56,139,253,0.5)")
        fig.add_annotation(x=L_ss/2, y=D*0.25, text="Liquid", font=dict(color="#388bfd"), showarrow=False)
        fig.add_annotation(x=L_ss/2, y=D*0.75, text="Gas", font=dict(color="#8b949e"), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0.05*L_ss, y0=0, x1=0.95*L_ss, y1=D*0.25,
                      line=dict(width=0), fillcolor="rgba(63,185,80,0.5)")
        fig.add_shape(type="rect", x0=0.05*L_ss, y0=D*0.25, x1=0.95*L_ss, y1=D*0.60,
                      line=dict(width=0), fillcolor="rgba(56,139,253,0.5)")
        fig.add_annotation(x=L_ss/2, y=D*0.12, text="Water", font=dict(color="#3fb950"), showarrow=False)
        fig.add_annotation(x=L_ss/2, y=D*0.42, text="Oil", font=dict(color="#388bfd"), showarrow=False)
        fig.add_annotation(x=L_ss/2, y=D*0.80, text="Gas", font=dict(color="#8b949e"), showarrow=False)

    # Dimension arrows
    fig.add_annotation(x=L_ss/2, y=-D*0.18,
                       text=f"L_ss = {L_ss:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=-0.03*L_ss, y=D/2,
                       text=f"D={D:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False, textangle=-90)

    fig.update_layout(
        **PLOTLY_THEME,
        title="Vessel Schematic (Side View)",
        height=260, margin=dict(l=60, r=20, t=40, b=40),
        showlegend=False,
    )
    fig.update_xaxes(visible=False, range=[-0.08*L_ss, 1.05*L_ss])
    fig.update_yaxes(visible=False, range=[-D*0.35, D*1.15], scaleanchor="x")
    return fig

def plot_vessel_schematic_vertical(D, H_total, H_liq, H_gas, n_phases=2, H_water=0, H_oil=0):
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, y0=0, x1=D, y1=H_total,
                  line=dict(color="#f0a500", width=2), fillcolor="#21262d")
    if n_phases == 2:
        fig.add_shape(type="rect", x0=0.05*D, y0=0, x1=0.95*D, y1=H_liq,
                      fillcolor="rgba(56,139,253,0.5)", line=dict(width=0))
        fig.add_annotation(x=D/2, y=H_liq/2, text="Liquid", font=dict(color="#388bfd"), showarrow=False)
        fig.add_annotation(x=D/2, y=H_liq + H_gas/2, text="Gas", font=dict(color="#8b949e"), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0.05*D, y0=0, x1=0.95*D, y1=H_water,
                      fillcolor="rgba(63,185,80,0.5)", line=dict(width=0))
        fig.add_shape(type="rect", x0=0.05*D, y0=H_water, x1=0.95*D, y1=H_water+H_oil,
                      fillcolor="rgba(56,139,253,0.5)", line=dict(width=0))
        fig.add_annotation(x=D/2, y=H_water/2, text="Water", font=dict(color="#3fb950"), showarrow=False)
        if H_oil > 0.05:
            fig.add_annotation(x=D/2, y=H_water+H_oil/2, text="Oil", font=dict(color="#388bfd"), showarrow=False)
        fig.add_annotation(x=D/2, y=H_water+H_oil+H_gas/2, text="Gas", font=dict(color="#8b949e"), showarrow=False)

    fig.add_annotation(x=D*1.2, y=H_total/2,
                       text=f"H={H_total:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=D/2, y=-H_total*0.06,
                       text=f"D={D:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)

    fig.update_layout(
        **PLOTLY_THEME,
        title="Vessel Schematic (Front View)",
        height=380, margin=dict(l=40, r=80, t=40, b=40),
        showlegend=False,
    )
    fig.update_xaxes(visible=False, range=[-0.1*D, 1.5*D])
    fig.update_yaxes(visible=False, range=[-H_total*0.12, H_total*1.1], scaleanchor="x")
    return fig



def plot_horizontal_with_nozzles(D, L_ss, nz_pos, weir, is_three):
    """Enhanced horizontal schematic showing nozzle positions and weir."""
    fig = go.Figure()

    # Shell
    fig.add_shape(type="rect", x0=0, y0=0, x1=L_ss, y1=D,
                  line=dict(color="#f0a500", width=2.5), fillcolor="#1a1f2e")

    # Phase zones
    if is_three:
        fig.add_shape(type="rect", x0=0, y0=0, x1=L_ss, y1=D*0.25,
                      line=dict(width=0), fillcolor="rgba(63,185,80,0.25)")
        fig.add_shape(type="rect", x0=0, y0=D*0.25, x1=L_ss, y1=D*0.60,
                      line=dict(width=0), fillcolor="rgba(56,139,253,0.20)")
        fig.add_annotation(x=L_ss*0.15, y=D*0.12, text="Water",
                           font=dict(color="#3fb950", size=10), showarrow=False)
        fig.add_annotation(x=L_ss*0.15, y=D*0.42, text="Oil",
                           font=dict(color="#388bfd", size=10), showarrow=False)
        fig.add_annotation(x=L_ss*0.15, y=D*0.80, text="Gas",
                           font=dict(color="#8b949e", size=10), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0, y0=0, x1=L_ss, y1=D*0.50,
                      line=dict(width=0), fillcolor="rgba(56,139,253,0.25)")
        fig.add_annotation(x=L_ss*0.15, y=D*0.25, text="Liquid",
                           font=dict(color="#388bfd", size=10), showarrow=False)
        fig.add_annotation(x=L_ss*0.15, y=D*0.75, text="Gas",
                           font=dict(color="#8b949e", size=10), showarrow=False)

    # Weir lines
    if is_three and weir:
        xw = weir["x_weir_m"]
        # Water weir
        fig.add_shape(type="line", x0=xw, y0=0, x1=xw, y1=weir["h_water_weir_m"],
                      line=dict(color="#3fb950", width=2.5, dash="solid"))
        fig.add_annotation(x=xw, y=weir["h_water_weir_m"]+D*0.03,
                           text=f"Water weir\n{weir['h_water_weir_m']:.3f} m",
                           font=dict(color="#3fb950", size=9), showarrow=False)
        # Oil overflow weir
        fig.add_shape(type="line", x0=xw+D*0.05, y0=0, x1=xw+D*0.05, y1=weir["h_oil_weir_m"],
                      line=dict(color="#388bfd", width=2.5, dash="solid"))
        fig.add_annotation(x=xw+D*0.05, y=weir["h_oil_weir_m"]+D*0.03,
                           text=f"Oil weir\n{weir['h_oil_weir_m']:.3f} m",
                           font=dict(color="#388bfd", size=9), showarrow=False)

    # Nozzle markers
    nozzle_colors = {
        "Inlet": "#f0a500", "Gas Outlet": "#8b949e", "Oil Outlet": "#388bfd",
        "Water Outlet": "#3fb950", "Liquid Outlet": "#388bfd",
        "Drain": "#ff6b6b", "Relief/Vent": "#da70d6", "Level Gauge": "#ffd700",
    }
    nozzle_symbols = {
        "Inlet": "arrow-right", "Gas Outlet": "arrow-up", "Oil Outlet": "arrow-right",
        "Water Outlet": "arrow-down", "Liquid Outlet": "arrow-down",
        "Drain": "triangle-down", "Relief/Vent": "arrow-up", "Level Gauge": "diamond",
    }
    for name, (px, py, desc) in nz_pos.items():
        col = nozzle_colors.get(name, "#ffffff")
        sym = nozzle_symbols.get(name, "circle")
        fig.add_trace(go.Scatter(
            x=[px], y=[py], mode="markers+text",
            marker=dict(symbol=sym, size=12, color=col,
                        line=dict(color="#0d1117", width=1)),
            text=[name], textposition="top center",
            textfont=dict(color=col, size=9),
            name=name, showlegend=True,
        ))
        # Leader line to vessel wall
        if py > D * 0.85:   # top nozzle — extend upward
            fig.add_shape(type="line", x0=px, y0=D, x1=px, y1=D+D*0.12,
                          line=dict(color=col, width=1, dash="dot"))
        elif py < D * 0.15:  # bottom nozzle — extend downward
            fig.add_shape(type="line", x0=px, y0=0, x1=px, y1=-D*0.12,
                          line=dict(color=col, width=1, dash="dot"))

    # Dimension labels
    fig.add_annotation(x=L_ss/2, y=-D*0.28, text=f"L_ss = {L_ss:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=-L_ss*0.04, y=D/2, text=f"D = {D:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False, textangle=-90)

    fig.update_layout(
        **PLOTLY_THEME,
        title="Horizontal Separator — Nozzle Layout",
        height=420,
        margin=dict(l=70, r=30, t=50, b=60),
        legend=dict(orientation="h", y=-0.25, font=dict(size=9),
                    bgcolor="rgba(0,0,0,0)"),
    )
    fig.update_xaxes(visible=False, range=[-L_ss*0.08, L_ss*1.06])
    fig.update_yaxes(visible=False, range=[-D*0.45, D*1.25], scaleanchor="x")
    return fig


def plot_vertical_with_nozzles(D, H_total, nz_pos, best, is_three):
    """Enhanced vertical schematic with nozzle positions."""
    fig = go.Figure()

    H_liq = best.get("H_liq_m", best.get("H_oil_m", 0) + best.get("H_water_m", 0))

    # Shell
    fig.add_shape(type="rect", x0=0, y0=0, x1=D, y1=H_total,
                  line=dict(color="#f0a500", width=2.5), fillcolor="#1a1f2e")

    # Phase zones
    if is_three:
        H_w = best.get("H_water_m", 0)
        H_o = best.get("H_oil_m", 0)
        fig.add_shape(type="rect", x0=0.05*D, y0=0, x1=0.95*D, y1=H_w,
                      fillcolor="rgba(63,185,80,0.30)", line=dict(width=0))
        fig.add_shape(type="rect", x0=0.05*D, y0=H_w, x1=0.95*D, y1=H_w+H_o,
                      fillcolor="rgba(56,139,253,0.25)", line=dict(width=0))
        fig.add_annotation(x=D/2, y=H_w/2, text="Water",
                           font=dict(color="#3fb950", size=9), showarrow=False)
        if H_o > 0.05:
            fig.add_annotation(x=D/2, y=H_w+H_o/2, text="Oil",
                               font=dict(color="#388bfd", size=9), showarrow=False)
        fig.add_annotation(x=D/2, y=H_w+H_o+best.get("H_gas_m",0.6)/2, text="Gas",
                           font=dict(color="#8b949e", size=9), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0.05*D, y0=0, x1=0.95*D, y1=H_liq,
                      fillcolor="rgba(56,139,253,0.25)", line=dict(width=0))
        fig.add_annotation(x=D/2, y=H_liq/2, text="Liquid",
                           font=dict(color="#388bfd", size=9), showarrow=False)

    # Demister pad
    H_mist_h = best.get("H_mist_m", me_thickness if "me_thickness" in dir() else 0.15)
    H_mist_bot = H_liq + best.get("H_gas_m", max(0.6, D*0.5))
    fig.add_shape(type="rect", x0=0.05*D, y0=H_mist_bot, x1=0.95*D, y1=H_mist_bot+H_mist_h,
                  fillcolor="rgba(240,165,0,0.35)", line=dict(color="#f0a500", width=1))
    fig.add_annotation(x=D*1.05, y=H_mist_bot+0.075, text="Demister",
                       font=dict(color="#f0a500", size=8), showarrow=False)

    # Nozzles
    nozzle_colors = {
        "Inlet": "#f0a500", "Gas Outlet": "#8b949e", "Oil Outlet": "#388bfd",
        "Water Outlet": "#3fb950", "Liquid Outlet": "#388bfd",
        "Drain": "#ff6b6b", "Relief/Vent": "#da70d6", "Level Gauge": "#ffd700",
    }
    for name, (px, pz, desc) in nz_pos.items():
        col = nozzle_colors.get(name, "#ffffff")
        if pz >= H_total * 0.9:   # top
            x_sym, xline0, xline1 = D/2, D/2, D/2
            y_sym = H_total
            fig.add_shape(type="line", x0=D/2, y0=H_total, x1=D/2, y1=H_total+D*0.15,
                          line=dict(color=col, width=1.5, dash="dot"))
        elif pz <= 0.01:   # bottom
            fig.add_shape(type="line", x0=D/2, y0=0, x1=D/2, y1=-D*0.15,
                          line=dict(color=col, width=1.5, dash="dot"))
            y_sym = 0
        else:   # side nozzle — on the right wall
            fig.add_shape(type="line", x0=D, y0=pz, x1=D+D*0.2, y1=pz,
                          line=dict(color=col, width=1.5, dash="dot"))
            y_sym = pz

        fig.add_trace(go.Scatter(
            x=[D*0.5], y=[y_sym if (pz <= 0.01 or pz >= H_total*0.9) else pz],
            mode="markers+text",
            marker=dict(symbol="circle", size=10, color=col,
                        line=dict(color="#0d1117", width=1)),
            text=[name], textposition="middle right",
            textfont=dict(color=col, size=9),
            name=name, showlegend=True,
        ))

    # Dimensions
    fig.add_annotation(x=D*1.35, y=H_total/2, text=f"H = {H_total:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=D/2, y=-D*0.25, text=f"D = {D:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)

    fig.update_layout(
        **PLOTLY_THEME,
        title="Vertical Separator — Nozzle Layout",
        height=520,
        margin=dict(l=40, r=100, t=50, b=60),
        legend=dict(orientation="v", x=1.05, font=dict(size=9),
                    bgcolor="rgba(0,0,0,0)"),
    )
    fig.update_xaxes(visible=False, range=[-D*0.3, D*1.8])
    fig.update_yaxes(visible=False, range=[-H_total*0.15, H_total*1.15], scaleanchor="x")
    return fig

# ─────────────────────────────────────────────────────────────────────────────
# SIDEBAR INPUTS
# ─────────────────────────────────────────────────────────────────────────────

with st.sidebar:
    # ── Logo ──────────────────────────────────────────────────────────────────
    col1, col2, col3 = st.columns([1,2,1])

    with col2:
        st.image("LOGO.png")
    st.markdown("---")

    st.markdown("### ⚡ Configuration")
    orientation = st.selectbox("Orientation", ["Horizontal", "Vertical"])
    n_phases = st.selectbox("Phases", ["Two-Phase (Gas/Liquid)", "Three-Phase (Gas/Oil/Water)"])
    is_three = "Three" in n_phases

    st.markdown("### 🌡 Operating Conditions")
    P_barg = st.number_input("Operating Pressure [barg]", 0.0, 200.0, 4.0, step=0.5,
                              help="Gauge pressure in bar. Converted to kPa abs internally (P_abs = (P_barg + 1.01325) × 100)")
    P_kPaa = (P_barg + 1.01325) * 100.0  # convert barg → kPa abs
    T_C = st.number_input("Operating Temperature [°C]", -50.0, 200.0, 15.0, step=1.0)
    T_K = T_C + 273.15

    st.markdown("### 💨 Gas Properties")
    MW_gas = st.number_input("Gas Molecular Weight [kg/kmol]", 10.0, 60.0, 22.0, step=0.5)
    Z_factor = st.number_input("Gas Compressibility Z [-]", 0.50, 1.20, 0.85, step=0.01)
    Qg_m3d = st.number_input("Gas Flow Rate [m³/d]", 1.0, 1e8, 35_000.0, step=10_000.0,
                              help="Standard cubic metres per day (0°C, 101.325 kPa)")
    # Convert to m³/s at operating conditions
    Qg_std = Qg_m3d / 86400   # m³/s standard
    P_std = 101.325  # kPa
    T_std = 273.15   # K
    rho_g_std = gas_density(P_std, T_std, MW_gas, 1.0)
    rho_g = gas_density(P_kPaa, T_K, MW_gas, Z_factor)
    mu_g = gas_viscosity_estimate(T_K, MW_gas)
    Qg_op = Qg_std * (P_std / P_kPaa) * (T_K / T_std) * Z_factor  # actual m³/s

    st.markdown("### 🛢 Oil / Liquid Properties")
    if is_three:
        API_grav = st.number_input("Oil API Gravity [°API]", 10.0, 60.0, 35.0, step=1.0)
        rho_o_15 = 141500 / (API_grav + 131.5)  # kg/m³ at 15°C
        rho_o = liquid_density_correction(rho_o_15, T_C)
        Qo_m3d = st.number_input("Oil Flow Rate [m³/d]", 0.1, 10000.0, 100.0, step=5.0)
        Qo = Qo_m3d / 86400
        mu_o = st.number_input("Oil Viscosity [mPa·s]", 0.5, 500.0, 8.0, step=0.5) * 1e-3

        rho_w = st.number_input("Water Density [kg/m³]", 990.0, 1100.0, 1050.0, step=5.0)
        Qw_m3d = st.number_input("Water Flow Rate [m³/d]", 0.1, 10000.0, 50.0, step=5.0)
        Qw = Qw_m3d / 86400
        mu_w = st.number_input("Water Viscosity [mPa·s]", 0.3, 3.0, 0.8, step=0.1) * 1e-3

        rho_l = rho_o  # representative liquid for gas section
        Ql = Qo + Qw
    else:
        API_grav = st.number_input("Liquid API Gravity [°API]", 10.0, 80.0, 40.0, step=1.0)
        rho_l_15 = 141500 / (API_grav + 131.5)
        rho_l = liquid_density_correction(rho_l_15, T_C)
        Ql_m3d = st.number_input("Liquid Flow Rate [m³/d]", 0.1, 10000.0, 80.0, step=5.0)
        Ql = Ql_m3d / 86400
        mu_o = None; mu_w = None; rho_o = rho_l; rho_w = None; Qo = Ql; Qw = 0

    st.markdown("### 🔬 Droplet Design")
    dp_mic = st.slider("Design Droplet Diameter [µm]", 50, 500, 150, step=10,
                        help="API 12J: 100–200 µm typical for mist extractors")
    if is_three:
        dp_w_mic = st.slider("Water Droplet in Oil [µm]", 50, 500, 200, step=10)
    else:
        dp_w_mic = 200

    st.markdown("### 🔩 Inlet Device")
    inlet_device = st.selectbox(
        "Inlet device type",
        list(INLET_DEVICES.keys()),
        index=0,
        help="Inlet device affects K factor, vessel length, and nozzle ρv² limit.",
    )
    dev_params = get_inlet_device_params(inlet_device)
    st.caption(dev_params["desc"])
    # Show effect summary as small metrics
    ic1, ic2, ic3 = st.columns(3)
    ic1.metric("K factor", f"×{dev_params['K_derating']:.2f}")
    ic2.metric("Length", f"×{dev_params['L_correction']:.2f}")
    ic3.metric("Pre-sep", f"{dev_params['pre_sep_credit']*100:.0f}%")

    st.markdown("### 📐 Souders-Brown K Factor")
    k_mode = st.radio(
        "K factor selection",
        ["Auto (pressure-based, GPSA Fig. 7-3)", "Preset (by internals type)", "Manual override"],
        label_visibility="collapsed",
    )
    if k_mode == "Auto (pressure-based, GPSA Fig. 7-3)":
        K_auto = souders_brown_K(T_K, P_kPaa)
        st.caption(f"Auto K = **{K_auto:.4f} m/s** at {P_barg:.1f} barg")
        K_user = None  # will use auto inside sizing functions
    elif k_mode == "Preset (by internals type)":
        preset = st.selectbox("Internals type", [
            "Wire mesh demister pad  (K = 0.107 m/s)",
            "Vane pack / chevron     (K = 0.122 m/s)",
            "No internals / gravity  (K = 0.061 m/s)",
            "High-capacity mesh pad  (K = 0.150 m/s)",
            "Agglomeration filter    (K = 0.046 m/s)",
        ])
        preset_map = {
            "Wire mesh demister pad  (K = 0.107 m/s)": 0.107,
            "Vane pack / chevron     (K = 0.122 m/s)": 0.122,
            "No internals / gravity  (K = 0.061 m/s)": 0.061,
            "High-capacity mesh pad  (K = 0.150 m/s)": 0.150,
            "Agglomeration filter    (K = 0.046 m/s)": 0.046,
        }
        K_user = preset_map[preset]
        st.caption(f"Preset K = **{K_user:.4f} m/s**")
    else:
        K_user = st.number_input("K factor [m/s]", 0.010, 0.300, 0.107, step=0.001,
                                  format="%.3f",
                                  help="Souders-Brown K [m/s]. Typical range: 0.046–0.150 m/s")
        st.caption(f"Manual K = **{K_user:.4f} m/s**")

    st.markdown("### 🌫️ Mist Extractor")
    me_type = st.selectbox(
        "Mist extractor type",
        list(MIST_EXTRACTORS.keys()),
        index=1,   # default: wire mesh pad (standard)
        help="Select mist extractor type. Affects K factor, vessel height, and outlet gas quality.",
    )
    me_data = MIST_EXTRACTORS[me_type]
    st.caption(me_data["desc"])
    if me_type != "None (gravity separation only)":
        mc1, mc2, mc3 = st.columns(3)
        mc1.metric("Cut size", f"{me_data['dp_micron']} µm")
        mc2.metric("Thickness", f"{me_data['thickness_m']*1000:.0f} mm")
        mc3.metric("K", f"{me_data['K_factor']} m/s")

    # Derive ME params for passing to sizing functions
    me_thickness = me_data["thickness_m"]
    me_K_derating = me_data["K_derating"]  # applied on top of inlet device derating

    st.markdown("### ⏱ Retention Times")
    if is_three:
        t_ret_o = retention_time_liquid(API_grav, 0, "three")
        t_ret_w = water_retention_time(T_C)
        tr_o_override = st.number_input("Oil Retention Time [min]", 1.0, 60.0,
                                         round(t_ret_o/60, 1), step=0.5)
        tr_w_override = st.number_input("Water Retention Time [min]", 1.0, 60.0,
                                         round(t_ret_w/60, 1), step=0.5)
        t_ret_o = tr_o_override * 60
        t_ret_w = tr_w_override * 60
        t_ret_l = t_ret_o
    else:
        t_ret_l = retention_time_liquid(API_grav, 0, "two")
        tr_override = st.number_input("Liquid Retention Time [min]", 1.0, 60.0,
                                       round(t_ret_l/60, 1), step=0.5)
        t_ret_l = tr_override * 60

    calc_btn = st.button("▶  CALCULATE", use_container_width=True, type="primary")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN LAYOUT
# ─────────────────────────────────────────────────────────────────────────────

st.markdown("""
<div style='display:flex;align-items:center;gap:1rem;margin-bottom:0.5rem'>
  <div>
    <h1 style='margin:0'>OpenSep | Separator Sizing Calculator</h1>
    <div class='section-pill'>API 12J  ·  GPSA EDB  ·  SI UNITS</div>
    <p style='color:#8b949e;font-size:0.85rem;margin:0'>
      Oil &amp; Gas Phase Separator — Two-Phase &amp; Three-Phase  |  Horizontal &amp; Vertical
    </p>
  </div>
</div>
<hr>
""", unsafe_allow_html=True)

tab_results, tab_table, tab_nozzles, tab_me, tab_theory, tab_params = st.tabs([
    "📊  Results", "📋  Sizing Table", "🔧  Nozzles & Weir", "🌫️  Mist Extractor", "📖  Theory & Equations", "🔢  Input Summary"
])

# ─────────────────────────────────────────────────────────────────────────────
if calc_btn or True:  # always display (updates on any sidebar change)

    # Run calculation
    best, df_all = None, None

    if orientation == "Horizontal" and not is_three:
        best, df_all = size_horizontal_two_phase(
            Qg_op, Ql, rho_g, rho_l, mu_g, dp_mic, t_ret_l, T_K, P_kPaa,
            K_override=K_user, inlet_device=inlet_device,
            me_K_derating=me_K_derating
        )
    elif orientation == "Horizontal" and is_three:
        best, df_all = size_horizontal_three_phase(
            Qg_op, Qo, Qw, rho_g, rho_o, rho_w,
            mu_g, mu_o, mu_w, dp_mic, dp_w_mic,
            t_ret_o, t_ret_w, T_K, P_kPaa,
            K_override=K_user, inlet_device=inlet_device,
            me_K_derating=me_K_derating
        )
    elif orientation == "Vertical" and not is_three:
        best, df_all = size_vertical_two_phase(
            Qg_op, Ql, rho_g, rho_l, mu_g, dp_mic, t_ret_l, T_K, P_kPaa,
            K_override=K_user, inlet_device=inlet_device,
            me_thickness=me_thickness, me_K_derating=me_K_derating
        )
    else:
        best, df_all = size_vertical_three_phase(
            Qg_op, Qo, Qw, rho_g, rho_o, rho_w,
            mu_g, mu_o, mu_w, dp_mic, dp_w_mic,
            t_ret_o, t_ret_w, T_K, P_kPaa,
            K_override=K_user, inlet_device=inlet_device,
            me_thickness=me_thickness, me_K_derating=me_K_derating
        )

    # ── TAB 1: RESULTS ────────────────────────────────────────────────────────
    with tab_results:
        if best is None:
            st.error("⚠️ No feasible vessel found. Check input flow rates and conditions.")
        else:
            D = best["D_m"]
            is_horiz = orientation == "Horizontal"

            # Primary metrics row
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.markdown(f"""
                <div class='result-card'>
                  <div class='label'>Vessel Diameter</div>
                  <div class='value'>{D:.3f}<span class='unit'>m</span></div>
                </div>""", unsafe_allow_html=True)
            with col2:
                if is_horiz:
                    st.markdown(f"""
                    <div class='result-card'>
                      <div class='label'>Seam-to-Seam Length</div>
                      <div class='value'>{best['L_ss_m']:.3f}<span class='unit'>m</span></div>
                    </div>""", unsafe_allow_html=True)
                else:
                    st.markdown(f"""
                    <div class='result-card'>
                      <div class='label'>Total Height</div>
                      <div class='value'>{best['H_total_m']:.3f}<span class='unit'>m</span></div>
                    </div>""", unsafe_allow_html=True)
            with col3:
                ld_val = best.get("LD", best.get("HD", 0))
                st.markdown(f"""
                <div class='result-card'>
                  <div class='label'>{'L/D' if is_horiz else 'H/D'} Ratio</div>
                  <div class='value'>{ld_val:.2f}</div>
                </div>""", unsafe_allow_html=True)
            with col4:
                vg = best["v_gas_ms"]
                vm = best["v_max_ms"]
                util = vg / vm * 100
                st.markdown(f"""
                <div class='result-card'>
                  <div class='label'>Gas Velocity Utilisation</div>
                  <div class='value'>{util:.1f}<span class='unit'>%</span></div>
                </div>""", unsafe_allow_html=True)

            st.markdown("<br>", unsafe_allow_html=True)

            # Status checks
            # ── ME effect diagnostic ─────────────────────────────────────
            # Determine which length is controlling
            if is_horiz and best is not None:
                _ctrl_keys = ["L_gas_ctrl_m","L_oil_ctrl_m","L_wat_ctrl_m"]
                _ctrl_vals = {k: best.get(k,0) for k in _ctrl_keys}
                _ctrl_vals["L_ret_o"] = best.get("V_oil_m3",0) / (0.35*math.pi*best["D_m"]**2/4) if is_three else 0
                _ctrl_vals["L_ret_w"] = best.get("V_wat_m3",0) / (0.25*math.pi*best["D_m"]**2/4) if is_three else 0
                _ctrl_vals["L_ret_liq"] = best.get("V_liq_m3",0) / (0.5*math.pi*best["D_m"]**2/4) if not is_three else 0
                _ctrl_name = max(_ctrl_vals, key=_ctrl_vals.get)
                _gas_controls = _ctrl_name == "L_gas_ctrl_m"
            else:
                _gas_controls = False

            if me_type == "None (gravity separation only)":
                if is_horiz and not _gas_controls:
                    st.markdown(
                        f'<div class="warn-card">⚠️ <strong>No mist extractor — vessel size unchanged.</strong><br>'
                        f'Controlling length is <strong>liquid-liquid separation</strong> (not gas settling), '
                        f'so the K factor derating has no effect on vessel dimensions. '
                        f'However: outlet gas quality is degraded (liquid carryover), '
                        f'the inlet nozzle ρv² limit remains 1490 kg/(m·s²), '
                        f'and this design is not suitable where gas outlet quality matters.</div>',
                        unsafe_allow_html=True)
                elif is_horiz and _gas_controls:
                    st.markdown(
                        '<div class="warn-card">⚠️ <strong>No mist extractor:</strong> '
                        'Gas settling controls — K derated 25%, vessel is larger than with ME.</div>',
                        unsafe_allow_html=True)
                else:
                    st.markdown(
                        '<div class="warn-card">⚠️ <strong>No mist extractor:</strong> '
                        'K derated 25%, gas disengagement zone increased to 0.9 m minimum. '
                        'Outlet gas quality degraded.</div>',
                        unsafe_allow_html=True)
            else:
                me_k = MIST_EXTRACTORS[me_type]["K_factor"]
                me_dp = MIST_EXTRACTORS[me_type]["dp_micron"]
                if is_horiz and not _gas_controls:
                    st.markdown(
                        f'<div class="ok-card">✅ <strong>Mist extractor: {me_type}</strong> — '
                        f'K={me_k} m/s, cut size={me_dp} µm. '
                        f'Note: vessel size is controlled by <strong>liquid retention / liquid-liquid settling</strong>, '
                        f'not gas capacity — so changing ME type does not change vessel dimensions in this case.</div>',
                        unsafe_allow_html=True)
                else:
                    st.markdown(
                        f'<div class="ok-card">✅ <strong>Mist extractor: {me_type}</strong> — '
                        f'K={me_k} m/s, cut size={me_dp} µm, '
                        f'thickness={MIST_EXTRACTORS[me_type]["thickness_m"]*1000:.0f} mm.</div>',
                        unsafe_allow_html=True)

            ld = best.get("LD", best.get("HD", 0))
            lo, hi = (3, 5) if is_horiz else (1.5, 4.0)
            if lo <= ld <= hi:
                st.markdown(f'<div class="ok-card">✅ L/D = {ld:.2f} — within API 12J preferred range ({lo}–{hi})</div>',
                            unsafe_allow_html=True)
            else:
                st.markdown(f'<div class="warn-card">⚠️ L/D = {ld:.2f} — outside preferred range ({lo}–{hi}). Consider different diameter.</div>',
                            unsafe_allow_html=True)

            if util > 85:
                st.markdown(f'<div class="warn-card">⚠️ Gas velocity at {util:.0f}% of Souders-Brown limit. Consider larger diameter.</div>',
                            unsafe_allow_html=True)
            else:
                st.markdown(f'<div class="ok-card">✅ Gas velocity at {util:.0f}% of Souders-Brown limit — acceptable.</div>',
                            unsafe_allow_html=True)

            st.markdown("<br>", unsafe_allow_html=True)
            lcol, rcol = st.columns([1.2, 1])

            with lcol:
                # Schematic
                if is_horiz:
                    fig_sch = plot_vessel_schematic_horizontal(
                        D, best["L_ss_m"], n_phases=3 if is_three else 2
                    )
                else:
                    hw = best.get("H_water_m", best.get("H_liq_m", 0))
                    ho = best.get("H_oil_m", best.get("H_liq_m", hw))
                    if not is_three:
                        ho = best.get("H_liq_m", 0)
                        hw = 0
                    fig_sch = plot_vessel_schematic_vertical(
                        D, best["H_total_m"],
                        H_liq=best.get("H_liq_m", ho + hw),
                        H_gas=best["H_gas_m"],
                        n_phases=3 if is_three else 2,
                        H_water=hw, H_oil=ho,
                    )
                st.plotly_chart(fig_sch, width='stretch')

            with rcol:
                # Detail results
                st.markdown("**Design Details**")
                details = {
                    "Gas Flow (actual)": f"{Qg_op*3600:.2f} m³/h",
                    "Gas Density": f"{rho_g:.3f} kg/m³",
                    "Gas Viscosity": f"{mu_g*1e6:.2f} µPa·s",
                    "Souders-Brown K": f"{best['K_factor']:.4f} m/s",
                    "Max Gas Velocity": f"{best['v_max_ms']:.4f} m/s",
                    "Actual Gas Velocity": f"{best['v_gas_ms']:.4f} m/s",
                }
                details["Inlet Device"] = best.get("inlet_device", inlet_device)
                details["K Derating"] = f"{best.get('K_derating', dev_params['K_derating']):.2f}×"
                details["Length Correction"] = f"{best.get('L_correction', dev_params['L_correction']):.2f}×"
                if best.get("pre_sep_credit", 0) > 0:
                    details["Pre-sep Credit"] = f"{best.get('pre_sep_credit', 0):.0f}%"
                if "vt_ms" in best:
                    details["Droplet Terminal Vel."] = f"{best['vt_ms']*1000:.3f} mm/s"
                if "efficiency" in best:
                    details["Droplet Removal Eff."] = f"{best['efficiency']:.1f} %"
                details["Mist Extractor"] = me_type
                if is_three and is_horiz:
                    details["Oil Ret. Volume"] = f"{best['V_oil_m3']:.3f} m³"
                    details["Water Ret. Volume"] = f"{best['V_wat_m3']:.3f} m³"
                else:
                    details["Liquid Ret. Volume"] = f"{best.get('V_liq_m3', Ql*t_ret_l):.3f} m³"

                for k, v in details.items():
                    c1, c2 = st.columns([2, 1.5])
                    c1.markdown(f"<span style='color:#8b949e;font-size:0.8rem'>{k}</span>",
                                unsafe_allow_html=True)
                    c2.markdown(f"<span style='color:#e6edf3;font-size:0.85rem;font-family:IBM Plex Mono'>{v}</span>",
                                unsafe_allow_html=True)

            # LD curve
            if df_all is not None and len(df_all) > 3:
                ld_col = "LD" if "LD" in df_all.columns else "HD"
                ylabel = "L/D Ratio" if is_horiz else "H/D Ratio"
                title = f"{'L/D' if is_horiz else 'H/D'} vs Vessel Diameter"
                fig_ld = plot_LD_curve(df_all, best, ylabel=ylabel, title=title)
                st.plotly_chart(fig_ld, width='stretch')

    # ── TAB 2: TABLE ─────────────────────────────────────────────────────────
    with tab_table:
        if df_all is not None:
            st.markdown('<div class="section-pill">All Feasible Vessel Sizes</div>', unsafe_allow_html=True)
            # Highlight the selected row
            def highlight_best(row):
                if abs(row["D_m"] - best["D_m"]) < 0.001:
                    return ["background-color: #1c1810; color: #f0a500"] * len(row)
                return [""] * len(row)
            st.dataframe(
                df_all.style.apply(highlight_best, axis=1).format(precision=3),
                width='stretch',
                height=480,
            )
            # Bar chart: L_eff vs D
            bar_y = "L_eff_m" if "L_eff_m" in df_all else "H_total_m"
            fig_bar = px.bar(df_all, x="D_m", y=bar_y,
                             color="D_m",
                             color_continuous_scale=[[0, "#21262d"], [1, "#f0a500"]],
                             labels={"D_m": "Diameter [m]", bar_y: "Length / Height [m]"},
                             title="Effective Length vs Diameter")
            fig_bar.add_vline(x=best["D_m"], line=dict(color="#f0a500", dash="dash"),
                              annotation_text="Selected")
            fig_bar.update_layout(**PLOTLY_THEME, height=300, margin=dict(l=50, r=20, t=40, b=40))
            st.plotly_chart(fig_bar, width='stretch')
        else:
            st.warning("No results to display.")


    # ── TAB 3: NOZZLES & WEIR ────────────────────────────────────────────────
    with tab_nozzles:
        if best is None:
            st.warning("Run a calculation first.")
        else:
            st.markdown('<div class="section-pill">API 12J §6 — Nozzle Sizing & Positions</div>', unsafe_allow_html=True)

            # ── Compute nozzle sizes ──────────────────────────────────────────
            # Inlet: combined gas + liquid
            Q_in_total = Qg_op + Ql
            rho_in = (Qg_op * rho_g + Ql * rho_l) / Q_in_total if Q_in_total > 0 else rho_g
            nz_inlet = nozzle_diameter(Q_in_total, rho_in, "inlet",
                                        rho_v2_limit=dev_params["rho_v2_limit"])

            # Gas outlet
            nz_gas = nozzle_diameter(Qg_op, rho_g, "gas_outlet")

            # Liquid outlet(s)
            if is_three:
                nz_oil = nozzle_diameter(Qo, rho_o, "liquid_outlet")
                nz_wat = nozzle_diameter(Qw, rho_w, "water_outlet")
            else:
                nz_liq = nozzle_diameter(Ql, rho_l, "liquid_outlet")

            # Drain (assume 5% of liquid flow for sizing)
            nz_drain = nozzle_diameter(Ql * 0.05, rho_l, "liquid_outlet")

            # ── Nozzle summary table ──────────────────────────────────────────
            # Inlet device summary banner
            st.markdown(f"""
<div class='ok-card' style='border-color:#f0a500;background:#1c1810;color:#f0a500'>
  <strong>🔩 Inlet Device: {inlet_device}</strong><br>
  <span style='font-size:0.8rem;color:#8b949e'>{dev_params['desc']}</span><br>
  K derating: <strong>{dev_params['K_derating']:.2f}×</strong> &nbsp;|&nbsp;
  Length correction: <strong>{dev_params['L_correction']:.2f}×</strong> &nbsp;|&nbsp;
  Pre-separation credit: <strong>{dev_params['pre_sep_credit']*100:.0f}%</strong> &nbsp;|&nbsp;
  Inlet ρv² limit: <strong>{dev_params['rho_v2_limit']} kg/(m·s²)</strong>
</div>""", unsafe_allow_html=True)

            st.markdown("#### Nozzle Schedule")

            nz_rows = []
            def dn_mm(d): return f"DN {int(round(d*1000))}"

            nz_rows.append({
                "Tag": "N1", "Service": "Inlet (Gas+Liquid)",
                "Flow [m³/s]": f"{Q_in_total:.4f}",
                "Fluid Density [kg/m³]": f"{rho_in:.1f}",
                "d_calc [mm]": f"{nz_inlet['d_calc_m']*1000:.1f}",
                "d_selected": dn_mm(nz_inlet['d_sel_m']),
                "Velocity [m/s]": f"{nz_inlet['velocity_ms']:.2f}",
                "Check Criterion": nz_inlet['criterion'],
            })
            nz_rows.append({
                "Tag": "N2", "Service": "Gas Outlet",
                "Flow [m³/s]": f"{Qg_op:.4f}",
                "Fluid Density [kg/m³]": f"{rho_g:.3f}",
                "d_calc [mm]": f"{nz_gas['d_calc_m']*1000:.1f}",
                "d_selected": dn_mm(nz_gas['d_sel_m']),
                "Velocity [m/s]": f"{nz_gas['velocity_ms']:.2f}",
                "Check Criterion": nz_gas['criterion'],
            })

            if is_three:
                nz_rows.append({
                    "Tag": "N3", "Service": "Oil Outlet",
                    "Flow [m³/s]": f"{Qo:.4f}",
                    "Fluid Density [kg/m³]": f"{rho_o:.1f}",
                    "d_calc [mm]": f"{nz_oil['d_calc_m']*1000:.1f}",
                    "d_selected": dn_mm(nz_oil['d_sel_m']),
                    "Velocity [m/s]": f"{nz_oil['velocity_ms']:.2f}",
                    "Check Criterion": nz_oil['criterion'],
                })
                nz_rows.append({
                    "Tag": "N4", "Service": "Water Outlet",
                    "Flow [m³/s]": f"{Qw:.4f}",
                    "Fluid Density [kg/m³]": f"{rho_w:.1f}",
                    "d_calc [mm]": f"{nz_wat['d_calc_m']*1000:.1f}",
                    "d_selected": dn_mm(nz_wat['d_sel_m']),
                    "Velocity [m/s]": f"{nz_wat['velocity_ms']:.2f}",
                    "Check Criterion": nz_wat['criterion'],
                })
            else:
                nz_rows.append({
                    "Tag": "N3", "Service": "Liquid Outlet",
                    "Flow [m³/s]": f"{Ql:.4f}",
                    "Fluid Density [kg/m³]": f"{rho_l:.1f}",
                    "d_calc [mm]": f"{nz_liq['d_calc_m']*1000:.1f}",
                    "d_selected": dn_mm(nz_liq['d_sel_m']),
                    "Velocity [m/s]": f"{nz_liq['velocity_ms']:.2f}",
                    "Check Criterion": nz_liq['criterion'],
                })

            nz_rows.append({
                "Tag": "N5", "Service": "Drain",
                "Flow [m³/s]": f"{Ql*0.05:.5f}",
                "Fluid Density [kg/m³]": f"{rho_l:.1f}",
                "d_calc [mm]": f"{nz_drain['d_calc_m']*1000:.1f}",
                "d_selected": dn_mm(nz_drain['d_sel_m']),
                "Velocity [m/s]": f"{nz_drain['velocity_ms']:.2f}",
                "Check Criterion": nz_drain['criterion'],
            })

            df_nozzles = pd.DataFrame(nz_rows)
            st.dataframe(df_nozzles, width="stretch", hide_index=True)

            st.markdown("---")

            # ── Nozzle positions ──────────────────────────────────────────────
            D_sel = best["D_m"]
            is_horiz_nz = orientation == "Horizontal"

            if is_horiz_nz:
                L_ss_sel = best["L_ss_m"]
                nz_pos = nozzle_positions_horizontal(D_sel, L_ss_sel, best, is_three)
                weir = weir_position_horizontal(
                    D_sel, L_ss_sel, best, is_three,
                    Qo_m3s=Qo if is_three else 0,
                    Qw_m3s=Qw if is_three else 0,
                    rho_o=rho_o, rho_w=rho_w if is_three else 1050
                )
            else:
                H_total_sel = best["H_total_m"]
                nz_pos = nozzle_positions_vertical(D_sel, H_total_sel, best, is_three)
                weir = None

            st.markdown("#### Nozzle Positions")
            pos_rows = []
            for name, (px, py, desc) in nz_pos.items():
                if is_horiz_nz:
                    pos_rows.append({"Nozzle": name, "x [m] (axial)": f"{px:.3f}",
                                     "y [m] (elevation)": f"{py:.3f}", "Notes": desc})
                else:
                    pos_rows.append({"Nozzle": name, "z [m] (elevation)": f"{py:.3f}", "Notes": desc})
            st.dataframe(pd.DataFrame(pos_rows), width="stretch", hide_index=True)

            # ── Weir (horizontal three-phase only) ───────────────────────────
            if is_horiz_nz and is_three and weir:
                st.markdown("---")
                st.markdown("#### Weir Positions & Heights")
                st.markdown('<div class="section-pill">API 12J §5.4.3 — Weir Design</div>', unsafe_allow_html=True)

                wc1, wc2, wc3 = st.columns(3)
                wc1.markdown(f"""
                <div class='result-card'>
                  <div class='label'>Weir Axial Position</div>
                  <div class='value'>{weir['x_weir_m']:.3f}<span class='unit'>m</span></div>
                </div>""", unsafe_allow_html=True)
                wc2.markdown(f"""
                <div class='result-card'>
                  <div class='label'>Water Weir Height</div>
                  <div class='value'>{weir['h_water_weir_m']:.3f}<span class='unit'>m</span></div>
                </div>""", unsafe_allow_html=True)
                wc3.markdown(f"""
                <div class='result-card'>
                  <div class='label'>Oil Overflow Weir Height</div>
                  <div class='value'>{weir['h_oil_weir_m']:.3f}<span class='unit'>m</span></div>
                </div>""", unsafe_allow_html=True)

                st.markdown(f"""
                <div class='ok-card'>
                  📐 {weir['desc']}<br>
                  Oil pad thickness between water weir and oil overflow weir:
                  <strong>{weir['oil_pad_thickness_m']:.3f} m</strong>
                </div>""", unsafe_allow_html=True)

            # ── Schematic with nozzles ────────────────────────────────────────
            st.markdown("---")
            st.markdown("#### Vessel Layout with Nozzles")

            if is_horiz_nz:
                fig_nz = plot_horizontal_with_nozzles(
                    D_sel, best["L_ss_m"], nz_pos, weir, is_three
                )
            else:
                fig_nz = plot_vertical_with_nozzles(
                    D_sel, best["H_total_m"], nz_pos, best, is_three
                )
            st.plotly_chart(fig_nz, width="stretch")

            # ── Basis note ────────────────────────────────────────────────────
            st.info(
                "**Nozzle sizing basis — API 12J §6:**\n"
                "- Inlet: ρv² ≤ 1490 kg/(m·s²)  \n"
                "- Gas outlet: v ≤ 18 m/s  \n"
                "- Liquid outlets: v ≤ 3 m/s  \n"
                "- Selected to next standard DN per ASME B36.10.  \n"
                "- Nozzle positions are indicative; final layout per vessel vendor drawing."
            )


    # ── TAB 4: MIST EXTRACTOR ────────────────────────────────────────────────
    with tab_me:
        if best is None:
            st.warning("Run a calculation first.")
        else:
            st.markdown('<div class="section-pill">API 12J §8 / GPSA Ch. 7 — Mist Extractor</div>',
                        unsafe_allow_html=True)

            D_sel = best["D_m"]
            K_sb = best.get("K_factor", souders_brown_K(T_K, P_kPaa))

            me_result = size_mist_extractor(
                me_type, Qg_op, rho_g, rho_l if not is_three else rho_o,
                D_sel, orientation=orientation, K_sb=K_sb
            )

            # ── No mist extractor ────────────────────────────────────────────
            if me_result["no_me"]:
                st.markdown(f"""
<div class='warn-card'>
  ⚠️ <strong>No mist extractor selected.</strong><br>
  {me_result['desc']}<br>
  K factor has been derated by <strong>25%</strong> in the separator sizing.
  Outlet gas will contain significant liquid carryover. Not recommended for
  sales gas, compression suction, or any duty with a gas quality requirement.
</div>""", unsafe_allow_html=True)

            else:
                # ── Status card ───────────────────────────────────────────────
                util = me_result["utilisation_pct"]
                if me_result.get("undersized"):
                    st.markdown(f"""
<div class='warn-card'>
  ⚠️ Mist extractor <strong>undersized</strong> for D={D_sel} m.
  Required area {me_result['A_required_m2']:.3f} m² > actual {me_result['A_actual_m2']:.3f} m².
  Consider increasing vessel diameter or switching to a higher-K extractor type.
</div>""", unsafe_allow_html=True)
                elif util > 85:
                    st.markdown(f"""
<div class='warn-card'>
  ⚠️ Gas velocity through extractor at {util:.1f}% of allowable — approaching flooding.
  Consider increasing vessel diameter.
</div>""", unsafe_allow_html=True)
                else:
                    st.markdown(f"""
<div class='ok-card'>
  ✅ Mist extractor sized OK. Gas velocity utilisation: <strong>{util:.1f}%</strong>
</div>""", unsafe_allow_html=True)

                # ── Primary metrics ───────────────────────────────────────────
                m1, m2, m3, m4 = st.columns(4)
                m1.markdown(f"""
<div class='result-card'>
  <div class='label'>Required Area</div>
  <div class='value'>{me_result['A_required_m2']:.3f}<span class='unit'>m²</span></div>
</div>""", unsafe_allow_html=True)
                m2.markdown(f"""
<div class='result-card'>
  <div class='label'>Actual Area</div>
  <div class='value'>{me_result['A_actual_m2']:.3f}<span class='unit'>m²</span></div>
</div>""", unsafe_allow_html=True)
                m3.markdown(f"""
<div class='result-card'>
  <div class='label'>Gas Velocity</div>
  <div class='value'>{me_result['v_me_ms']:.3f}<span class='unit'>m/s</span></div>
</div>""", unsafe_allow_html=True)
                m4.markdown(f"""
<div class='result-card'>
  <div class='label'>Utilisation</div>
  <div class='value'>{me_result['utilisation_pct']:.1f}<span class='unit'>%</span></div>
</div>""", unsafe_allow_html=True)

                st.markdown("<br>", unsafe_allow_html=True)

                # ── Properties table ──────────────────────────────────────────
                st.markdown("#### Mist Extractor Properties")
                prop_rows = [
                    ("Type", me_result["me_type"]),
                    ("Design K factor", f"{me_result['K_design']:.4f} m/s"),
                    ("Allowable gas velocity", f"{me_result['v_allow_ms']:.4f} m/s"),
                    ("Droplet cut size (d₅₀)", f"{me_result['dp_cut_micron']} µm" if me_result['dp_cut_micron'] else "—"),
                    ("Pad / element thickness", f"{me_result['thickness_m']*1000:.0f} mm"),
                    ("Wire diameter", f"{me_result['wire_dia_mm']} mm" if me_result['wire_dia_mm'] else "—"),
                    ("Void fraction", f"{me_result['void_fraction']:.2f}" if me_result['void_fraction'] else "—"),
                    ("Specific surface area", f"{me_result['specific_area']} m²/m³" if me_result['specific_area'] else "—"),
                    ("ΔP clean", f"{me_result['dP_clean_Pa']:.0f} Pa  ({me_result['dP_clean_Pa']/1000:.3f} kPa)"),
                    ("ΔP loaded (design)", f"{me_result['dP_loaded_Pa']:.0f} Pa  ({me_result['dP_loaded_Pa']/1000:.3f} kPa)"),
                    ("Orientation", orientation),
                    ("Vessel diameter", f"{D_sel:.3f} m"),
                ]
                df_me = pd.DataFrame(prop_rows, columns=["Parameter", "Value"])
                st.dataframe(df_me, width="stretch", hide_index=True)

                # ── Pressure drop chart ───────────────────────────────────────
                st.markdown("#### Pressure Drop vs Gas Velocity")
                v_range = np.linspace(0.01, me_result["v_allow_ms"] * 1.2, 80)
                # Simplified dP ∝ v² (Souders-Brown correlation style)
                dP_range = me_result["dP_clean_Pa"] * (v_range / me_result["v_me_ms"])**2
                fig_dp = go.Figure()
                fig_dp.add_trace(go.Scatter(
                    x=v_range, y=dP_range, mode="lines",
                    line=dict(color="#388bfd", width=2), name="Est. ΔP"
                ))
                fig_dp.add_vline(x=me_result["v_me_ms"],
                                  line=dict(color="#f0a500", dash="dash", width=1.5),
                                  annotation_text=f"Operating {me_result['v_me_ms']:.3f} m/s",
                                  annotation_font=dict(color="#f0a500", size=10))
                fig_dp.add_vline(x=me_result["v_allow_ms"],
                                  line=dict(color="#ff6b6b", dash="dot", width=1.5),
                                  annotation_text=f"Max {me_result['v_allow_ms']:.3f} m/s",
                                  annotation_font=dict(color="#ff6b6b", size=10))
                fig_dp.add_hrect(y0=me_result["dP_clean_Pa"], y1=me_result["dP_loaded_Pa"],
                                  fillcolor="rgba(240,165,0,0.08)",
                                  line=dict(color="rgba(240,165,0,0.3)"),
                                  annotation_text="Clean→Loaded range",
                                  annotation_font=dict(color="#f0a500", size=9))
                fig_dp.update_layout(
                    **PLOTLY_THEME, **_AXIS_DEFAULTS,
                    title="Mist Extractor ΔP vs Gas Velocity (estimated)",
                    xaxis_title="Gas Velocity [m/s]",
                    yaxis_title="Pressure Drop [Pa]",
                    height=300, margin=dict(l=60, r=20, t=45, b=40),
                )
                st.plotly_chart(fig_dp, width="stretch")

                # ── Position note ──────────────────────────────────────────────
                st.markdown("#### Installation Notes")
                if orientation == "Vertical":
                    H_mist_bot = (best.get("H_liq_m",
                                  best.get("H_oil_m",0)+best.get("H_water_m",0))
                                  + best.get("H_gas_m", 0.6))
                    H_mist_top = H_mist_bot + me_result["thickness_m"]
                    st.markdown(f"""
- **Position:** {H_mist_bot:.3f} m to {H_mist_top:.3f} m above bottom tangent line
- Mounted horizontally across full vessel cross-section
- Support ring to leave ≥ 50 mm annular clearance at vessel wall
- Minimum 150 mm clearance below gas outlet nozzle
""")
                else:
                    st.markdown(f"""
- Mounted in **gas space** (upper {100*0.45:.0f}% of cross-section)
- Available area: {me_result['A_actual_m2']:.3f} m² (45% of vessel cross-section)
- For horizontal vessels, mesh pad spans full length of gas zone
- Minimum 150 mm clearance from gas outlet nozzle
""")

            # ── Comparison table ──────────────────────────────────────────────
            st.markdown("---")
            st.markdown("#### Mist Extractor Type Comparison")
            cmp_rows = []
            for name, props in MIST_EXTRACTORS.items():
                if name == "None (gravity separation only)":
                    cmp_rows.append({
                        "Type": name,
                        "K [m/s]": "—",
                        "Cut size [µm]": "—",
                        "Thickness [mm]": "—",
                        "ΔP clean [Pa]": "0",
                        "ΔP loaded [Pa]": "0",
                        "Best for": "Low-spec, onshore, no outlet quality req."
                    })
                else:
                    cmp_rows.append({
                        "Type": name,
                        "K [m/s]": str(props["K_factor"]),
                        "Cut size [µm]": str(props["dp_micron"]),
                        "Thickness [mm]": str(int(props["thickness_m"]*1000)),
                        "ΔP clean [Pa]": str(int(props["dP_clean_Pa"])),
                        "ΔP loaded [Pa]": str(int(props["dP_loaded_Pa"])),
                        "Best for": props["desc"].split(".")[0],
                    })
            st.dataframe(pd.DataFrame(cmp_rows), width="stretch", hide_index=True)

    # ── TAB 5: THEORY ────────────────────────────────────────────────────────
    with tab_theory:
        st.markdown("## Theory & Reference Equations")
        st.markdown("""
### 1. Fluid Properties

**Gas Density** (Real Gas Law):
$$\\rho_g = \\frac{P \\cdot M_w}{Z \\cdot R \\cdot T}$$
where $P$ [Pa], $M_w$ [kg/kmol], $Z$ [-], $R = 8314$ J/(kmol·K), $T$ [K].

---

### 2. Droplet Settling Velocity

**Stokes' Law** (Re < 0.4):
$$v_t = \\frac{d_p^2 (\\rho_L - \\rho_G) g}{18 \\mu_G}$$

**General Drag Correlation** (Schiller-Naumann, API 12J / GPSA Ch. 7):
$$v_t = \\sqrt{\\frac{4 d_p (\\rho_L - \\rho_G) g}{3 C_D \\rho_G}}$$

where the drag coefficient:
- $Re < 0.4$: $C_D = 24/Re$
- $0.4 \\le Re \\le 500$: $C_D = \\dfrac{24}{Re} + \\dfrac{6}{1+\\sqrt{Re}} + 0.4$
- $Re > 500$: $C_D = 0.44$

---

### 3. Souders-Brown Equation (API 12J Table 1 / GPSA Fig. 7-3)

$$v_{SB} = K \\sqrt{\\frac{\\rho_L - \\rho_G}{\\rho_G}}$$

The $K$ factor [m/s] is read from the API 12J / GPSA pressure-dependent chart.  
For separators with demister pads, $K$ ≈ 0.061–0.122 m/s depending on operating pressure.

---

### 4. Horizontal Two-Phase Sizing (API 12J §5.3)

**Minimum effective length from gas settling:**
$$L_{eff} \\geq \\frac{v_G \\cdot D}{v_t}$$

**Minimum effective length from liquid retention:**
$$L_{eff} \\geq \\frac{Q_L \\cdot t_R}{A_L}$$

where $A_L = f_L \\cdot \\pi D^2/4$ (liquid cross-section fraction, typically 50%).

**Seam-to-seam length:**
$$L_{ss} = L_{eff} + 0.3D \\quad \\text{(head allowance)}$$

---

### 5. Vertical Two-Phase Sizing (API 12J §5.2)

**Required gas area:**
$$A_{req} = \\frac{Q_G}{v_{design}}, \\quad v_{design} = \\min(v_{SB},\\; 0.75\\, v_t)$$

**Total height:**
$$H_T = H_{liq} + H_{gas} + H_{mist} + H_{inlet}$$

with $H_{gas} \\geq 0.6$ m (API 12J minimum), $H_{mist} = 0.15$ m, $H_{inlet} = 0.3$ m.

---

### 6. Three-Phase Sizing (API 12J §5.4)

**Liquid–liquid separation** — oil/water interface:

Oil droplet rise velocity in water:
$$v_{t,o/w} = \\frac{d_{p,o}^2 (\\rho_W - \\rho_O) g}{18 \\mu_W}$$

Water droplet settling in oil:
$$v_{t,w/o} = \\frac{d_{p,w}^2 (\\rho_W - \\rho_O) g}{18 \\mu_O}$$

**Retention times** (API 12J Table 2):

| Oil API Gravity | Oil Retention | Water Retention |
|----------------|--------------|----------------|
| > 35 °API      | 3 min        | 5–20 min       |
| 25–35 °API     | 5 min        | 10–20 min      |
| < 25 °API      | 10 min       | 20+ min        |

---

### 7. Design Criteria (API 12J / GPSA)

| Parameter | Horizontal | Vertical |
|-----------|-----------|---------|
| Preferred L/D or H/D | 3 – 5 | 1.5 – 4 |
| Max gas velocity | ≤ $v_{SB}$ | ≤ $\\min(v_{SB}, 0.75v_t)$ |
| Design droplet | 100–200 µm | 100–200 µm |
| Demister K factor | GPSA Fig. 7-3 | GPSA Fig. 7-3 |

---

### 8. Mist Extractor Sizing (API 12J §8 / GPSA Ch. 7)

The mist extractor removes entrained liquid droplets from the gas before it exits.

**Required cross-sectional area** (Souders-Brown):
$$A_{ME} = \\frac{Q_G}{K_{ME} \\sqrt{\\frac{\\rho_L - \\rho_G}{\\rho_G}}}$$

**Gas velocity through extractor:**
$$v_{ME} = \\frac{Q_G}{A_{ME,actual}}$$

**Utilisation** (must be < 100%, recommended < 85%):
$$U = \\frac{v_{ME}}{v_{allow}} \\times 100\\%$$

**Pressure drop** estimated as proportional to $v^2$:
$$\\Delta P = \\Delta P_{clean} \\left(\\frac{v_{ME}}{v_{ref}}\\right)^2$$

**Without mist extractor:** K factor derated by 25% per GPSA recommendation.

| Type | K [m/s] | d₅₀ [µm] | ΔP clean [Pa] | Notes |
|---|---|---|---|---|
| Wire mesh (std) | 0.107 | 10 | 125 | API 12J baseline |
| Wire mesh (hi-cap) | 0.122 | 20 | 100 | High gas rate |
| Vane pack | 0.122 | 30 | 200 | Slug tolerant |
| Axial cyclone | 0.150 | 5 | 500 | Offshore / high-spec |
| Agglomeration | 0.046 | 3 | 250 | Critical duty |

---

### 9. Inlet Device Selection (API 12J §6.3 / GPSA Ch. 7)

The inlet device controls how the incoming multiphase stream enters the vessel.
It affects three key parameters:

| Device | ρv² limit [kg/(m·s²)] | K derating | Length correction | Pre-sep credit |
|---|---|---|---|---|
| Half-pipe / Diverter plate | 1 490 | ×1.00 | ×1.00 | 0% |
| Perforated pipe / Spreader | 2 500 | ×0.95 | ×0.95 | 10% |
| Vane-type (Schoepentoeter) | 6 000 | ×1.05 | ×0.85 | 30% |
| Cyclonic inlet device | 12 000 | ×1.10 | ×0.75 | 60% |
| No device (nozzle only) | 1 490 | ×0.80 | ×1.30 | 0% |

**Pre-separation credit** reduces the effective liquid load entering the settling section:
$$Q_{eff} = Q_{inlet} \\times (1 - f_{pre-sep})$$

**Length correction** is applied to $L_{eff}$:
$$L_{eff,corrected} = L_{eff,calc} \\times f_{L}$$

**K derating** modifies the Souders-Brown allowable velocity:
$$K_{eff} = K_{base} \\times f_K$$

---

### 9. Nozzle Sizing (API 12J §6)

**Inlet nozzle** — momentum criterion (API 12J Eq. 6-1):
$$\\rho v^2 \\leq 1490 \\; \\text{kg/(m·s}^2\\text{)}$$
$$A_{inlet} \\geq Q_{total} \\sqrt{\\frac{\\rho}{1490}} \\quad \\Rightarrow \\quad d = \\sqrt{\\frac{4A}{\\pi}}$$

**Gas outlet**: $v \\leq 18$ m/s  
**Liquid outlet**: $v \\leq 3$ m/s  
**Relief/Vent**: $v \\leq 18$ m/s

Selected to next standard DN size per **ASME B36.10**.

---

### 9. Weir Design — Horizontal Three-Phase (API 12J §5.4.3)

The **water weir** controls the water/oil interface level:
$$h_{water\\,weir} = f_{water} \\cdot D \\quad (\\text{typically 25\\% of } D)$$

The **oil overflow weir** controls the oil level and prevents carry-over:
$$h_{oil\\,weir} = f_{oil} \\cdot D \\quad (\\text{typically 60\\% of } D)$$

**Oil pad thickness**:
$$\\Delta h_{oil} = h_{oil\\,weir} - h_{water\\,weir}$$

The weir is positioned at $x_{weir} = L_{eff}$ (end of the liquid retention zone).

---

### References
- **API 12J** — Specification for Oil and Gas Separators, 8th Edition
- **GPSA Engineering Data Book** — Chapter 7: Gas/Liquid Separators
- **API RP 12L** — Specification for Vertical and Horizontal Emulsion Treaters
- Schiller & Naumann (1935) drag correlation
- Souders & Brown (1934) demisting velocity
""")

    # ── TAB 4: INPUT SUMMARY ─────────────────────────────────────────────────
    with tab_params:
        st.markdown("## Input Parameter Summary")
        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown("**Operating Conditions**")
            params_op = {
                "Pressure": f"{P_barg:.2f} barg  ({P_kPaa:.1f} kPa abs)",
                "Temperature": f"{T_C:.1f} °C  ({T_K:.2f} K)",
                "Orientation": orientation,
                "Phase Configuration": "3-Phase" if is_three else "2-Phase",
            }
            for k, v in params_op.items():
                st.markdown(f"`{k}`: **{v}**")

            st.markdown("**Gas**")
            params_gas = {
                "Mol. Weight": f"{MW_gas:.1f} kg/kmol",
                "Z factor": f"{Z_factor:.3f}",
                "Density (op.)": f"{rho_g:.4f} kg/m³",
                "Viscosity (est.)": f"{mu_g*1e6:.2f} µPa·s",
                "Flow (std)": f"{Qg_m3d:,.0f} m³/d",
                "Flow (op.)": f"{Qg_op*3600:.4f} m³/h",
            }
            for k, v in params_gas.items():
                st.markdown(f"`{k}`: **{v}**")

        with col_b:
            if is_three:
                st.markdown("**Oil**")
                params_oil = {
                    "API Gravity": f"{API_grav:.1f} °API",
                    "Density (op.)": f"{rho_o:.2f} kg/m³",
                    "Viscosity": f"{mu_o*1000:.2f} mPa·s",
                    "Flow": f"{Qo*86400:.1f} m³/d",
                    "Retention Time": f"{t_ret_o/60:.1f} min",
                }
                for k, v in params_oil.items():
                    st.markdown(f"`{k}`: **{v}**")
                st.markdown("**Water**")
                params_w = {
                    "Density": f"{rho_w:.1f} kg/m³",
                    "Viscosity": f"{mu_w*1000:.2f} mPa·s",
                    "Flow": f"{Qw*86400:.1f} m³/d",
                    "Retention Time": f"{t_ret_w/60:.1f} min",
                }
                for k, v in params_w.items():
                    st.markdown(f"`{k}`: **{v}**")
            else:
                st.markdown("**Liquid**")
                params_liq = {
                    "API Gravity": f"{API_grav:.1f} °API",
                    "Density (op.)": f"{rho_l:.2f} kg/m³",
                    "Flow": f"{Ql*86400:.1f} m³/d",
                    "Retention Time": f"{t_ret_l/60:.1f} min",
                }
                for k, v in params_liq.items():
                    st.markdown(f"`{k}`: **{v}**")

            st.markdown("**Droplet Design**")
            st.markdown(f"`Gas section droplet`: **{dp_mic} µm**")
            if is_three:
                st.markdown(f"`Water-in-oil droplet`: **{dp_w_mic} µm**")

        st.info("ℹ️ All calculations are performed in SI units (Pa, K, m, kg, s). "
                "Results are for preliminary sizing only. Final design must comply with "
                "applicable codes (ASME, PED, API 12J) and be reviewed by a qualified engineer.")

# ─────────────────────────────────────────────────────────────────────────────
# FOOTER
# ─────────────────────────────────────────────────────────────────────────────
st.markdown("""
<hr>
<p style='text-align:center;color:#8b949e;font-size:0.72rem;font-family:IBM Plex Mono'>
API 12J · GPSA Engineering Data Book · SI Units Throughout · For Preliminary Sizing Only
</p>
""", unsafe_allow_html=True)