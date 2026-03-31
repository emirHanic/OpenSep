"""
calculations.py — OpenSep
All engineering calculation functions.
No Streamlit imports. Pure Python + math + numpy.
Follows API 12J (8th Ed.) and GPSA Engineering Data Book.
"""

import math
import numpy as np
import pandas as pd

from constants import (
    INLET_DEVICES, MIST_EXTRACTORS, STANDARD_NOZZLE_DN,
    HORIZ_3P_FRACTIONS, HORIZ_2P_LIQ_FRAC,
    R_GAS, G_ACCEL, P_STD, T_STD,
)

# ─────────────────────────────────────────────────────────────────────────────
# FLUID PROPERTIES
# ─────────────────────────────────────────────────────────────────────────────

def gas_density(P_kPa, T_K, MW, Z=1.0):
    """ρ_g [kg/m³] from real gas law.  R = 8314 J/(kmol·K)."""
    return (P_kPa * 1000.0 * MW) / (Z * R_GAS * T_K)


def liquid_density_correction(rho_l_15, T_C):
    """Thermal expansion correction for crude (API 11.1 light version).
    Δρ ≈ -0.65 kg/m³ per °C above 15 °C."""
    return rho_l_15 - 0.65 * (T_C - 15.0)


def gas_viscosity_estimate(T_K, MW):
    """Estimate gas viscosity [Pa·s] — Chapman-Enskog simplified."""
    mu = (2.6693e-6 * math.sqrt(MW * T_K)) / (3.5 ** 2 * 1.22)
    return max(mu, 1e-6)


def cunningham_correction(dp_m, P_kPa, T_K, MW_gas):
    """Cunningham slip correction factor (relevant for dp < 3 µm)."""
    mu  = gas_viscosity_estimate(T_K, MW_gas)
    rho = gas_density(P_kPa, T_K, MW_gas)
    lam = (mu / (0.499 * rho)) * math.sqrt(math.pi * MW_gas / (2.0 * R_GAS * T_K))
    Kn  = 2.0 * lam / dp_m
    return 1.0 + Kn * (1.257 + 0.4 * math.exp(-1.1 / Kn))


def actual_gas_flow(Qg_std_m3s, P_kPa, T_K, Z):
    """Convert standard gas flow [m³/s] to actual operating conditions."""
    return Qg_std_m3s * (P_STD / P_kPa) * (T_K / T_STD) * Z


# ─────────────────────────────────────────────────────────────────────────────
# SETTLING / TERMINAL VELOCITY
# ─────────────────────────────────────────────────────────────────────────────

def terminal_velocity_stokes(dp_m, rho_disp, rho_cont, mu_cont):
    """Stokes' law terminal velocity [m/s]."""
    return (dp_m**2 * abs(rho_disp - rho_cont) * G_ACCEL) / (18.0 * mu_cont)


def terminal_velocity_cd(dp_m, rho_disp, rho_cont, mu_cont):
    """
    Iterative Schiller-Naumann drag coefficient terminal velocity [m/s].
    Works for gas/liquid AND liquid/liquid systems.
    API 12J §5 / GPSA Ch. 7.
    """
    delta_rho = abs(rho_disp - rho_cont)
    if delta_rho < 1e-3 or mu_cont <= 0 or dp_m <= 0:
        return 1e-6

    vt = (dp_m**2 * delta_rho * G_ACCEL) / (18.0 * mu_cont)  # Stokes start
    for _ in range(200):
        Re = max(rho_cont * vt * dp_m / mu_cont, 1e-10)
        if Re < 0.4:
            CD = 24.0 / Re
        elif Re < 500.0:
            CD = 24.0 / Re + 6.0 / (1.0 + math.sqrt(Re)) + 0.4
        else:
            CD = 0.44
        vt_new = math.sqrt((4.0 * delta_rho * G_ACCEL * dp_m) / (3.0 * CD * rho_cont))
        if abs(vt_new - vt) / (vt + 1e-12) < 1e-6:
            return vt_new
        vt = vt_new
    return vt


# ─────────────────────────────────────────────────────────────────────────────
# SOUDERS-BROWN
# ─────────────────────────────────────────────────────────────────────────────

def souders_brown_K(T_K, P_kPa, service="gas_liquid"):
    """
    Pressure-dependent K [m/s] — GPSA Fig. 7-3 / API 12J Table 1.
    """
    P_psia = P_kPa * 0.145038
    if   P_psia <= 100:  K = 0.122
    elif P_psia <= 200:  K = 0.116
    elif P_psia <= 400:  K = 0.107
    elif P_psia <= 600:  K = 0.098
    elif P_psia <= 900:  K = 0.085
    elif P_psia <= 1200: K = 0.073
    else:                K = 0.061
    return K


def max_gas_velocity(K, rho_g, rho_l):
    """Allowable gas velocity [m/s] from Souders-Brown equation."""
    return K * math.sqrt((rho_l - rho_g) / rho_g)


# ─────────────────────────────────────────────────────────────────────────────
# RETENTION TIMES  (API 12J Table 2)
# ─────────────────────────────────────────────────────────────────────────────

def retention_time_liquid(API_gravity, GOR=0, phase="two"):
    """Recommended liquid retention time [s] per API 12J Table 2."""
    if   API_gravity > 35: t_r = 3 * 60
    elif API_gravity > 25: t_r = 5 * 60
    else:                  t_r = 10 * 60
    if phase == "three":
        t_r = max(t_r, 10 * 60)
    return t_r


def water_retention_time(T_C):
    """Water retention time [s] per API 12J."""
    if   T_C > 60: return  5 * 60
    elif T_C > 40: return 10 * 60
    else:          return 20 * 60


# ─────────────────────────────────────────────────────────────────────────────
# INLET DEVICE
# ─────────────────────────────────────────────────────────────────────────────

def get_inlet_device_params(device_name):
    """Return inlet device parameter dict."""
    return INLET_DEVICES.get(device_name, INLET_DEVICES["Half-pipe / Diverter plate"])


# ─────────────────────────────────────────────────────────────────────────────
# MIST EXTRACTOR SIZING  (API 12J §8 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

def size_mist_extractor(me_type, Qg_m3s, rho_g, rho_l, D_vessel,
                        orientation="Horizontal", K_sb=None):
    """
    Size the mist extractor per API 12J §8 / GPSA Ch. 7.
    Returns a sizing dict with area, velocity, utilisation, and dP.
    """
    me = MIST_EXTRACTORS.get(me_type, MIST_EXTRACTORS["Wire mesh pad (standard)"])

    if me_type == "None (gravity separation only)":
        return {
            "me_type": me_type, "desc": me["desc"],
            "K_design": None, "v_allow_ms": None,
            "A_required_m2": None, "A_actual_m2": None,
            "v_me_ms": None, "utilisation_pct": None,
            "dP_clean_Pa": 0.0, "dP_loaded_Pa": 0.0,
            "thickness_m": 0.0, "dp_cut_micron": None,
            "wire_dia_mm": None, "void_fraction": None,
            "specific_area": None, "K_derating": me["K_derating"],
            "undersized": False, "no_me": True,
        }

    K_design = me["K_factor"]
    if K_sb is not None:
        K_design = min(K_design, K_sb)

    v_allow    = K_design * math.sqrt((rho_l - rho_g) / rho_g)
    A_required = Qg_m3s / v_allow
    A_vessel   = math.pi * D_vessel**2 / 4
    A_actual   = 0.90 * A_vessel if orientation == "Vertical" else 0.45 * A_vessel
    v_me       = Qg_m3s / A_actual
    utilisation = (v_me / v_allow) * 100.0

    return {
        "me_type": me_type, "desc": me["desc"],
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
        "undersized": A_actual < A_required,
        "no_me": False,
    }


# ─────────────────────────────────────────────────────────────────────────────
# NOZZLE SIZING  (API 12J §6)
# ─────────────────────────────────────────────────────────────────────────────

def next_standard_nozzle(d_calc):
    """Return next standard DN nozzle size [m] >= d_calc."""
    for d in STANDARD_NOZZLE_DN:
        if d >= d_calc:
            return d
    return STANDARD_NOZZLE_DN[-1]


def nozzle_diameter(Q_m3s, rho, nozzle_type="liquid", rho_v2_limit=1490):
    """
    Size a nozzle per API 12J §6 velocity / momentum criteria.
    Returns dict with d_calc, d_selected, velocity, criterion.
    """
    if Q_m3s <= 0 or rho <= 0:
        return {"d_calc_m": 0.0, "d_sel_m": STANDARD_NOZZLE_DN[0],
                "velocity_ms": 0.0, "momentum": 0.0, "criterion": "—"}

    v_limits = {"gas_outlet": 18.0, "liquid_outlet": 3.0,
                "water_outlet": 3.0, "relief": 18.0}

    if nozzle_type == "inlet":
        rho_v2_limit = float(rho_v2_limit)
        A_min  = Q_m3s * math.sqrt(rho / rho_v2_limit)
        d_calc = math.sqrt(4.0 * A_min / math.pi)
        d_sel  = next_standard_nozzle(d_calc)
        v      = Q_m3s / (math.pi * d_sel**2 / 4)
        return {"d_calc_m": d_calc, "d_sel_m": d_sel, "velocity_ms": v,
                "momentum": rho * v**2,
                "criterion": f"ρv²={rho*v**2:.0f} ≤ {rho_v2_limit:.0f} kg/(m·s²)"}

    v_lim  = v_limits.get(nozzle_type, 3.0)
    d_calc = math.sqrt(4.0 * Q_m3s / (v_lim * math.pi))
    d_sel  = next_standard_nozzle(d_calc)
    v      = Q_m3s / (math.pi * d_sel**2 / 4)
    return {"d_calc_m": d_calc, "d_sel_m": d_sel, "velocity_ms": v,
            "momentum": rho * v**2, "criterion": f"v={v:.2f} ≤ {v_lim} m/s"}


# ─────────────────────────────────────────────────────────────────────────────
# NOZZLE & WEIR POSITIONS
# ─────────────────────────────────────────────────────────────────────────────

def nozzle_positions_horizontal(D, L_ss, best, is_three):
    """
    Nozzle centreline positions for a horizontal separator.
    Returns {name: (x_axial, y_elevation, description)}.
    x=0 at inlet head seam, y=0 at vessel bottom.
    """
    head_h    = 0.15 * D
    x_in      = head_h
    x_out     = L_ss - head_h

    pos = {
        "Inlet":       (x_in,  D * 0.65, f"x={x_in:.3f} m, y={D*0.65:.3f} m"),
        "Gas Outlet":  (x_out, D * 0.92, f"x={x_out:.3f} m, y={D*0.92:.3f} m (top)"),
        "Drain":       (L_ss * 0.5, 0.0, f"x={L_ss*0.5:.3f} m, bottom"),
        "Relief/Vent": (x_in + 0.3, D * 0.95, f"x={x_in+0.3:.3f} m, top"),
        "Level Gauge": (x_in + 0.2, D * 0.50, f"x={x_in+0.2:.3f} m, NLL"),
    }
    if is_three:
        pos["Oil Outlet"]   = (x_out, D * 0.60, f"x={x_out:.3f} m, y={D*0.60:.3f} m")
        pos["Water Outlet"] = (x_out, D * 0.08, f"x={x_out:.3f} m, y={D*0.08:.3f} m")
    else:
        pos["Liquid Outlet"] = (x_out, D * 0.15, f"x={x_out:.3f} m, y={D*0.15:.3f} m")
    return pos


def weir_position_horizontal(D, L_ss, best, is_three, Qo_m3s=0, Qw_m3s=0,
                              rho_o=800, rho_w=1050):
    """
    Weir heights and axial position for horizontal three-phase separator.
    API 12J §5.4.3.
    """
    if not is_three:
        return None
    x_weir      = best.get("L_eff_m", L_ss * 0.85)
    h_water_weir = D * 0.25
    h_oil_weir   = D * 0.60
    return {
        "x_weir_m":           round(x_weir, 3),
        "h_water_weir_m":     round(h_water_weir, 3),
        "h_oil_weir_m":       round(h_oil_weir, 3),
        "oil_pad_thickness_m": round(h_oil_weir - h_water_weir, 3),
        "desc": (
            f"Water weir at x={x_weir:.3f} m, h={h_water_weir:.3f} m.  "
            f"Oil overflow weir h={h_oil_weir:.3f} m.  "
            f"Oil pad={h_oil_weir - h_water_weir:.3f} m."
        ),
    }


def nozzle_positions_vertical(D, H_total, best, is_three):
    """
    Nozzle positions for a vertical separator.
    Returns {name: (x, z_elevation, description)}.
    z=0 at bottom tangent line.
    """
    H_liq       = best.get("H_liq_m", best.get("H_oil_m", 0) + best.get("H_water_m", 0))
    H_gas       = best.get("H_gas_m", max(0.6, D * 0.5))
    z_inlet     = H_liq + 0.15
    z_gas_out   = H_total - 0.1
    z_nll       = H_liq * 0.5

    pos = {
        "Inlet":       (D * 0.5, z_inlet,   f"z={z_inlet:.3f} m, side nozzle"),
        "Gas Outlet":  (D * 0.5, z_gas_out, f"z={z_gas_out:.3f} m (top)"),
        "Drain":       (D * 0.5, 0.0,        "z=0.0 m, bottom"),
        "Relief/Vent": (D * 0.5, z_gas_out,  f"z={z_gas_out:.3f} m, top"),
        "Level Gauge": (D * 0.5, z_nll,      f"z={z_nll:.3f} m, NLL"),
    }
    if is_three:
        z_oil   = best.get("H_water_m", 0) + best.get("H_oil_m", 0) * 0.5
        z_water = best.get("H_water_m", 0) * 0.15
        pos["Oil Outlet"]   = (D * 0.5, z_oil,   f"z={z_oil:.3f} m")
        pos["Water Outlet"] = (D * 0.5, z_water, f"z={z_water:.3f} m")
    else:
        pos["Liquid Outlet"] = (D * 0.5, H_liq * 0.1, f"z={H_liq*0.1:.3f} m")
    return pos


# ─────────────────────────────────────────────────────────────────────────────
# HORIZONTAL SEPARATOR SIZING  (API 12J §5.3 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

def _build_K(K_override, T_K, P_kPa, dev, me_K_derating):
    """Compose effective K from override / auto × inlet device × ME derating."""
    K_base = K_override if K_override is not None else souders_brown_K(T_K, P_kPa)
    return K_base * dev["K_derating"] * me_K_derating


def size_horizontal_two_phase(
    Qg_m3s, Ql_m3s,
    rho_g, rho_l, mu_g,
    dp_mic, t_ret_s,
    T_K, P_kPa,
    K_override=None,
    inlet_device="Half-pipe / Diverter plate",
    me_K_derating=1.0,
):
    """
    Horizontal two-phase separator sizing.
    API 12J §5.3 / GPSA Ch. 7.
    Returns (best_dict, full_dataframe).
    """
    dp_m  = dp_mic * 1e-6
    dev   = get_inlet_device_params(inlet_device)
    vt    = terminal_velocity_cd(dp_m, rho_l, rho_g, mu_g)
    K     = _build_K(K_override, T_K, P_kPa, dev, me_K_derating)
    v_max = max_gas_velocity(K, rho_g, rho_l)
    V_liq = Ql_m3s * t_ret_s
    f_liq = HORIZ_2P_LIQ_FRAC

    results = []
    for D in np.arange(0.3, 5.05, 0.05):
        A      = math.pi * D**2 / 4
        A_gas  = (1 - f_liq) * A
        A_liq  = f_liq * A
        v_gas  = Qg_m3s / A_gas
        if v_gas > v_max:
            continue

        L_gas_factor = 1.0 if me_K_derating >= 1.0 else 1.0 / me_K_derating
        L_min_gas    = (v_gas * D / vt) * L_gas_factor
        L_min_liq    = V_liq / A_liq
        L_eff        = max(L_min_gas, L_min_liq) * dev["L_correction"]
        LD           = L_eff / D
        if LD < 1.5:
            continue

        L_ss       = L_eff + 0.3 * D
        efficiency = min(1.0, vt * L_eff / (v_gas * D))
        results.append({
            "D_m": round(D, 2), "L_eff_m": round(L_eff, 3),
            "L_ss_m": round(L_ss, 3), "LD": round(LD, 2),
            "v_gas_ms": round(v_gas, 3), "v_max_ms": round(v_max, 3),
            "vt_ms": round(vt, 5), "efficiency": round(efficiency * 100, 1),
            "V_liq_m3": round(V_liq, 3), "K_factor": round(K, 4),
            "inlet_device": inlet_device, "K_derating": round(dev["K_derating"], 3),
            "L_correction": round(dev["L_correction"], 3),
        })

    return _select_best(results, ratio_key="LD", lo=3.0, hi=5.0)


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
    API 12J §5.4 / GPSA Ch. 7.
    """
    dp_o     = dp_o_mic * 1e-6
    dp_w     = dp_w_mic * 1e-6
    dev      = get_inlet_device_params(inlet_device)
    pre_sep  = dev["pre_sep_credit"]

    vt_gas   = terminal_velocity_cd(dp_o, rho_o, rho_g, mu_g)
    vt_o_w   = terminal_velocity_cd(dp_o, rho_o, rho_w, mu_w)
    vt_w_o   = terminal_velocity_cd(dp_w, rho_w, rho_o, mu_o)

    K        = _build_K(K_override, T_K, P_kPa, dev, me_K_derating)
    v_max    = max_gas_velocity(K, rho_g, rho_o)

    Qo_eff   = Qo_m3s * (1.0 - pre_sep)
    Qw_eff   = Qw_m3s * (1.0 - pre_sep)
    V_o      = Qo_eff * t_ret_o_s
    V_w      = Qw_eff * t_ret_w_s

    f        = HORIZ_3P_FRACTIONS
    results  = []
    for D in np.arange(0.3, 5.55, 0.05):
        A     = math.pi * D**2 / 4
        A_gas = f["gas"] * A
        A_oil = f["oil"] * A
        A_wat = f["water"] * A

        v_gas = Qg_m3s / A_gas
        if v_gas > v_max:
            continue

        L_gas_f = 1.0 if me_K_derating >= 1.0 else 1.0 / me_K_derating
        L_gas   = (v_gas * D / vt_gas)
        L_oilw  = (Qo_eff / A_oil) * D / vt_w_o
        L_wato  = (Qw_eff / A_wat) * D / vt_o_w
        L_ret_o = V_o / A_oil
        L_ret_w = V_w / A_wat
        L_eff   = max(L_gas * L_gas_f, L_oilw, L_wato, L_ret_o, L_ret_w) * dev["L_correction"]
        LD      = L_eff / D
        if LD < 1.5:
            continue

        L_ss = L_eff + 0.3 * D
        results.append({
            "D_m": round(D, 2), "L_eff_m": round(L_eff, 3),
            "L_ss_m": round(L_ss, 3), "LD": round(LD, 2),
            "v_gas_ms": round(v_gas, 4), "v_max_ms": round(v_max, 3),
            "L_gas_ctrl_m": round(L_gas, 3), "L_oil_ctrl_m": round(L_oilw, 3),
            "L_wat_ctrl_m": round(L_wato, 3),
            "V_oil_m3": round(V_o, 3), "V_wat_m3": round(V_w, 3),
            "K_factor": round(K, 4), "vt_ms": round(vt_gas, 5),
            "inlet_device": inlet_device, "K_derating": round(dev["K_derating"], 3),
            "L_correction": round(dev["L_correction"], 3),
            "pre_sep_credit": round(pre_sep * 100, 0),
        })

    return _select_best(results, ratio_key="LD", lo=3.0, hi=5.0)


# ─────────────────────────────────────────────────────────────────────────────
# VERTICAL SEPARATOR SIZING  (API 12J §5.2 / GPSA Ch. 7)
# ─────────────────────────────────────────────────────────────────────────────

def size_vertical_two_phase(
    Qg_m3s, Ql_m3s,
    rho_g, rho_l, mu_g,
    dp_mic, t_ret_s,
    T_K, P_kPa,
    K_override=None,
    inlet_device="Half-pipe / Diverter plate",
    me_thickness=0.15,
    me_K_derating=1.0,
):
    """Vertical two-phase separator sizing.  API 12J §5.2."""
    dp_m  = dp_mic * 1e-6
    dev   = get_inlet_device_params(inlet_device)
    vt    = terminal_velocity_cd(dp_m, rho_l, rho_g, mu_g)
    K     = _build_K(K_override, T_K, P_kPa, dev, me_K_derating)
    v_max = max_gas_velocity(K, rho_g, rho_l)
    v_des = min(v_max, 0.75 * vt)
    V_liq = Ql_m3s * t_ret_s

    results = []
    for D in np.arange(0.3, 5.05, 0.05):
        A     = math.pi * D**2 / 4
        v_gas = Qg_m3s / A
        if v_gas > v_max:
            continue
        H_liq     = V_liq / A
        H_gas_min = 0.6 if me_thickness > 0 else 0.9
        H_gas     = max(H_gas_min, D * 0.5)
        H_mist    = me_thickness
        H_total   = H_liq + H_gas + H_mist + 0.3   # +0.3 inlet zone
        HD        = H_total / D
        results.append({
            "D_m": round(D, 2), "H_liq_m": round(H_liq, 3),
            "H_gas_m": round(H_gas, 3), "H_mist_m": round(H_mist, 3),
            "H_total_m": round(H_total, 3), "HD": round(HD, 2),
            "v_gas_ms": round(v_gas, 4), "v_max_ms": round(v_max, 4),
            "vt_ms": round(vt, 5), "V_liq_m3": round(V_liq, 3),
            "K_factor": round(K, 4), "inlet_device": inlet_device,
            "K_derating": round(dev["K_derating"], 3),
        })

    return _select_best(results, ratio_key="HD", lo=1.5, hi=4.0)


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
    """Vertical three-phase separator sizing.  API 12J §5.4."""
    dp_o    = dp_o_mic * 1e-6
    dp_w    = dp_w_mic * 1e-6
    dev     = get_inlet_device_params(inlet_device)
    pre_sep = dev["pre_sep_credit"]

    vt_gas  = terminal_velocity_cd(dp_o, rho_o, rho_g, mu_g)
    vt_o_w  = terminal_velocity_cd(dp_o, rho_o, rho_w, mu_w)

    K       = _build_K(K_override, T_K, P_kPa, dev, me_K_derating)
    v_max   = max_gas_velocity(K, rho_o, rho_g) if rho_g < rho_o else K * math.sqrt((rho_o - rho_g) / rho_g)
    # Correct Souders-Brown for vertical
    v_max   = K * math.sqrt((rho_o - rho_g) / rho_g)

    V_o     = Qo_m3s * (1.0 - pre_sep) * t_ret_o_s
    V_w     = Qw_m3s * (1.0 - pre_sep) * t_ret_w_s

    results = []
    for D in np.arange(0.3, 5.55, 0.05):
        A     = math.pi * D**2 / 4
        v_gas = Qg_m3s / A
        if v_gas > v_max:
            continue
        H_oil     = V_o / A
        H_water   = V_w / A
        H_gas_min = 0.6 if me_thickness > 0 else 0.9
        H_gas     = max(H_gas_min, D * 0.5)
        H_mist    = me_thickness
        H_total   = H_water + H_oil + H_gas + H_mist + 0.3
        HD        = H_total / D
        results.append({
            "D_m": round(D, 2), "H_oil_m": round(H_oil, 3),
            "H_water_m": round(H_water, 3), "H_gas_m": round(H_gas, 3),
            "H_mist_m": round(H_mist, 3), "H_total_m": round(H_total, 3),
            "HD": round(HD, 2), "v_gas_ms": round(v_gas, 4),
            "v_max_ms": round(v_max, 4), "vt_ms": round(vt_gas, 5),
            "vt_gas_ms": round(vt_gas, 5), "vt_o_w_ms": round(vt_o_w, 5),
            "V_oil_m3": round(V_o, 3), "V_wat_m3": round(V_w, 3),
            "K_factor": round(K, 4),
        })

    return _select_best(results, ratio_key="HD", lo=1.5, hi=4.0)


# ─────────────────────────────────────────────────────────────────────────────
# INTERNAL HELPER
# ─────────────────────────────────────────────────────────────────────────────

def _select_best(results, ratio_key, lo, hi):
    """Pick preferred vessel: smallest satisfying lo ≤ ratio ≤ hi."""
    if not results:
        return None, None
    preferred = [r for r in results if lo <= r[ratio_key] <= hi]
    if not preferred:
        preferred = [r for r in results if r[ratio_key] <= hi + 1.0]
    if not preferred:
        preferred = results
    return preferred[0], pd.DataFrame(results)
