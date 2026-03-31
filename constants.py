"""
constants.py — OpenSep
All lookup tables and data constants.
Add new equipment types, standards, or data here without touching calculations.
"""

# ─────────────────────────────────────────────────────────────────────────────
# NOZZLE STANDARDS
# ─────────────────────────────────────────────────────────────────────────────

# Standard nozzle OD sizes [m] per ASME B36.10 / API 6D (DN 50 → DN 600)
STANDARD_NOZZLE_DN = [
    0.0508, 0.0635, 0.0762, 0.1016, 0.1524, 0.2032,
    0.2540, 0.3048, 0.3556, 0.4064, 0.5080, 0.6096,
]

# ─────────────────────────────────────────────────────────────────────────────
# INLET DEVICES  (API 12J §6.3 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────
INLET_DEVICES = {
    "Half-pipe / Diverter plate": {
        "rho_v2_limit":   1490,   # kg/(m·s²)  API 12J Eq. 6-1
        "K_derating":     1.00,   # no derating — baseline device
        "L_correction":   1.00,   # no length penalty
        "pre_sep_credit": 0.0,    # fraction of liquid pre-separated
        "desc": "Standard baseline device. API 12J ρv² ≤ 1490 kg/(m·s²).",
    },
    "Vane-type inlet (e.g. Schoepentoeter)": {
        "rho_v2_limit":   6000,
        "K_derating":     1.05,
        "L_correction":   0.85,
        "pre_sep_credit": 0.30,
        "desc": "High-capacity vane inlet. Higher ρv² limit, 15% L reduction, 30% pre-sep credit.",
    },
    "Cyclonic inlet device": {
        "rho_v2_limit":   12000,
        "K_derating":     1.10,
        "L_correction":   0.75,
        "pre_sep_credit": 0.60,
        "desc": "Cyclonic inlet. Highest capacity, 25% L reduction, 60% liquid pre-separation.",
    },
    "Inlet nozzle only (no device)": {
        "rho_v2_limit":   1490,
        "K_derating":     0.80,
        "L_correction":   1.30,
        "pre_sep_credit": 0.0,
        "desc": "No inlet device. K derated 20%, vessel length increased 30% per GPSA.",
    },
    "Perforated pipe / Spreader": {
        "rho_v2_limit":   2500,
        "K_derating":     0.95,
        "L_correction":   0.95,
        "pre_sep_credit": 0.10,
        "desc": "Perforated pipe distributor. Modest improvements over bare nozzle.",
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# MIST EXTRACTORS  (API 12J §8 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────
MIST_EXTRACTORS = {
    "None (gravity separation only)": {
        "K_factor":           None,
        "K_derating":         0.75,   # 25% K derating without mist extractor
        "dp_micron":          None,
        "thickness_m":        0.0,
        "dP_clean_Pa":        0.0,
        "dP_loaded_Pa":       0.0,
        "wire_diameter_mm":   None,
        "void_fraction":      None,
        "specific_area_m2m3": None,
        "desc": "No mist extractor. K derated 25%. Suitable only where outlet gas quality is not critical.",
    },
    "Wire mesh pad (standard)": {
        "K_factor":           0.107,
        "K_derating":         1.00,
        "dp_micron":          10.0,
        "thickness_m":        0.15,
        "dP_clean_Pa":        125.0,
        "dP_loaded_Pa":       375.0,
        "wire_diameter_mm":   0.28,
        "void_fraction":      0.97,
        "specific_area_m2m3": 200,
        "desc": "Standard knitted wire mesh. dp_cut ≈ 10 µm. GPSA K=0.107 m/s. API 12J §8.3.",
    },
    "Wire mesh pad (high capacity)": {
        "K_factor":           0.122,
        "K_derating":         1.00,
        "dp_micron":          20.0,
        "thickness_m":        0.20,
        "dP_clean_Pa":        100.0,
        "dP_loaded_Pa":       300.0,
        "wire_diameter_mm":   0.38,
        "void_fraction":      0.98,
        "specific_area_m2m3": 130,
        "desc": "High-capacity coarser mesh. Higher K, lower dP, larger cut size. For high gas rates.",
    },
    "Vane pack (chevron)": {
        "K_factor":           0.122,
        "K_derating":         1.00,
        "dp_micron":          30.0,
        "thickness_m":        0.30,
        "dP_clean_Pa":        200.0,
        "dP_loaded_Pa":       500.0,
        "wire_diameter_mm":   None,
        "void_fraction":      0.90,
        "specific_area_m2m3": 250,
        "desc": "Vane/chevron pack. Higher dP, handles slugging better than mesh. dp_cut ≈ 30 µm.",
    },
    "Axial cyclone bundle": {
        "K_factor":           0.150,
        "K_derating":         1.00,
        "dp_micron":          5.0,
        "thickness_m":        0.50,
        "dP_clean_Pa":        500.0,
        "dP_loaded_Pa":       1200.0,
        "wire_diameter_mm":   None,
        "void_fraction":      None,
        "specific_area_m2m3": None,
        "desc": "Axial cyclone tubes. Best efficiency (dp_cut ≈ 5 µm), highest K and dP. Offshore / high-spec.",
    },
    "Agglomeration filter": {
        "K_factor":           0.046,
        "K_derating":         1.00,
        "dp_micron":          3.0,
        "thickness_m":        0.25,
        "dP_clean_Pa":        250.0,
        "dP_loaded_Pa":       700.0,
        "wire_diameter_mm":   None,
        "void_fraction":      0.92,
        "specific_area_m2m3": 500,
        "desc": "Agglomeration / coalescing filter. Sub-micron capability. Lowest K; use for critical duty.",
    },
}

# ─────────────────────────────────────────────────────────────────────────────
# VESSEL GEOMETRY FRACTIONS  (API 12J Fig 5-2)
# ─────────────────────────────────────────────────────────────────────────────
# Cross-section fractions for horizontal three-phase
HORIZ_3P_FRACTIONS = {"gas": 0.40, "oil": 0.35, "water": 0.25}
# Liquid fill fraction for horizontal two-phase
HORIZ_2P_LIQ_FRAC = 0.50

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────
R_GAS    = 8314.0   # J/(kmol·K)
G_ACCEL  = 9.81     # m/s²
P_STD    = 101.325  # kPa  (standard conditions)
T_STD    = 273.15   # K
