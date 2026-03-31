"""
app.py — OpenSep
Streamlit UI entry point.
All engineering logic lives in calculations.py / constants.py.
All plot functions live in plots.py.
All CSS / HTML strings live in styles.py.
"""

import math
import streamlit as st
import numpy as np
import pandas as pd

# ── Local modules ─────────────────────────────────────────────────────────────
from constants import INLET_DEVICES, MIST_EXTRACTORS
from calculations import (
    gas_density, gas_viscosity_estimate, liquid_density_correction,
    actual_gas_flow, souders_brown_K,
    retention_time_liquid, water_retention_time,
    get_inlet_device_params, size_mist_extractor,
    nozzle_diameter, nozzle_positions_horizontal, nozzle_positions_vertical,
    weir_position_horizontal,
    size_horizontal_two_phase, size_horizontal_three_phase,
    size_vertical_two_phase,   size_vertical_three_phase,
)
from plots import (
    plot_ld_curve,
    plot_schematic_horizontal, plot_schematic_vertical,
    plot_horizontal_with_nozzles, plot_vertical_with_nozzles,
    plot_me_dp_curve, plot_length_vs_diameter,
)
from styles import CSS, HEADER_HTML, FOOTER_HTML, THEORY_MD
from styles import result_card, warn_card, ok_card, section_pill

# ─────────────────────────────────────────────────────────────────────────────
# PAGE CONFIG
# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="OpenSep",
    page_icon="⚙️",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.markdown(CSS, unsafe_allow_html=True)

# ─────────────────────────────────────────────────────────────────────────────
# SIDEBAR
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
st.markdown(HEADER_HTML, unsafe_allow_html=True)

tab_results, tab_table, tab_nozzles, tab_me, tab_theory, tab_params = st.tabs([
    "📊  Results", "📋  Sizing Table", "🔧  Nozzles & Weir",
    "🌫️  Mist Extractor", "📖  Theory & Equations", "🔢  Input Summary",
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
                    fig_sch = plot_schematic_horizontal(
                        D, best["L_ss_m"], n_phases=3 if is_three else 2
                    )
                else:
                    hw = best.get("H_water_m", best.get("H_liq_m", 0))
                    ho = best.get("H_oil_m", best.get("H_liq_m", hw))
                    if not is_three:
                        ho = best.get("H_liq_m", 0)
                        hw = 0
                    fig_sch = plot_schematic_vertical(
                        D, best["H_total_m"],
                        H_liq=best.get("H_liq_m", ho + hw),
                        H_gas=best["H_gas_m"],
                        n_phases=3 if is_three else 2,
                        H_water=hw, H_oil=ho,
                        H_mist=best.get('H_mist_m', me_thickness),
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
                fig_ld = plot_ld_curve(df_all, best, ylabel=ylabel, title=title)
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
            fig_bar = plot_length_vs_diameter(df_all, best)
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
                    D_sel, best["H_total_m"], nz_pos, best, is_three,
                    me_thickness=me_thickness
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
                fig_dp = plot_me_dp_curve(me_result)
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
        st.markdown(THEORY_MD)

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
st.markdown(FOOTER_HTML, unsafe_allow_html=True)
