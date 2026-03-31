"""
plots.py — OpenSep
All Plotly figure functions.
No Streamlit imports. Returns go.Figure objects.
"""

import math
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

# ─────────────────────────────────────────────────────────────────────────────
# THEME
# ─────────────────────────────────────────────────────────────────────────────
THEME = dict(
    paper_bgcolor="#0d1117",
    plot_bgcolor="#161b22",
    font=dict(family="IBM Plex Mono", color="#8b949e", size=11),
    title_font=dict(color="#e6edf3", size=13),
)
GRID = dict(
    xaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
    yaxis=dict(gridcolor="#21262d", linecolor="#30363d"),
)
AMBER  = "#f0a500"
BLUE   = "#388bfd"
GREEN  = "#3fb950"
GREY   = "#8b949e"
RED    = "#ff6b6b"
PURPLE = "#da70d6"
GOLD   = "#ffd700"

NOZZLE_COLORS = {
    "Inlet": AMBER, "Gas Outlet": GREY, "Oil Outlet": BLUE,
    "Water Outlet": GREEN, "Liquid Outlet": BLUE,
    "Drain": RED, "Relief/Vent": PURPLE, "Level Gauge": GOLD,
}
NOZZLE_SYMBOLS_H = {
    "Inlet": "arrow-right", "Gas Outlet": "arrow-up",
    "Oil Outlet": "arrow-right", "Water Outlet": "arrow-down",
    "Liquid Outlet": "arrow-down", "Drain": "triangle-down",
    "Relief/Vent": "arrow-up", "Level Gauge": "diamond",
}


# ─────────────────────────────────────────────────────────────────────────────
# L/D or H/D ENVELOPE CURVE
# ─────────────────────────────────────────────────────────────────────────────

def plot_ld_curve(df, best, ylabel="L/D Ratio", title="L/D vs Diameter"):
    ratio_col = "LD" if "LD" in df.columns else "HD"
    lo, hi = (3, 5) if ratio_col == "LD" else (1.5, 4)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["D_m"], y=df[ratio_col], mode="lines",
        line=dict(color=BLUE, width=2), name=ratio_col,
    ))
    fig.add_hrect(y0=lo, y1=hi, fillcolor="rgba(240,165,0,0.08)",
                  line=dict(color="rgba(240,165,0,0.25)"),
                  annotation_text="API 12J preferred",
                  annotation_font=dict(color=AMBER, size=10))
    fig.add_vline(x=best["D_m"], line=dict(color=AMBER, dash="dash", width=1.5),
                  annotation_text=f"D={best['D_m']} m",
                  annotation_font=dict(color=AMBER, size=10))
    fig.update_layout(**THEME, **GRID, title=title,
                      xaxis_title="D [m]", yaxis_title=ylabel,
                      height=320, margin=dict(l=50, r=20, t=40, b=40))
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# SIMPLE VESSEL SCHEMATICS (Results tab)
# ─────────────────────────────────────────────────────────────────────────────

def plot_schematic_horizontal(D, L_ss, n_phases=2):
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, y0=0, x1=L_ss, y1=D,
                  line=dict(color=AMBER, width=2), fillcolor="#21262d")
    if n_phases == 2:
        fig.add_shape(type="rect", x0=0.05*L_ss, y0=0, x1=0.95*L_ss, y1=D*0.50,
                      line=dict(width=0), fillcolor="rgba(56,139,253,0.5)")
        fig.add_annotation(x=L_ss/2, y=D*0.25, text="Liquid",
                           font=dict(color=BLUE), showarrow=False)
        fig.add_annotation(x=L_ss/2, y=D*0.75, text="Gas",
                           font=dict(color=GREY), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0.05*L_ss, y0=0, x1=0.95*L_ss, y1=D*0.25,
                      line=dict(width=0), fillcolor="rgba(63,185,80,0.5)")
        fig.add_shape(type="rect", x0=0.05*L_ss, y0=D*0.25, x1=0.95*L_ss, y1=D*0.60,
                      line=dict(width=0), fillcolor="rgba(56,139,253,0.5)")
        for y, text, col in [(D*0.12, "Water", GREEN), (D*0.42, "Oil", BLUE),
                              (D*0.80, "Gas", GREY)]:
            fig.add_annotation(x=L_ss/2, y=y, text=text,
                               font=dict(color=col), showarrow=False)
    fig.add_annotation(x=L_ss/2, y=-D*0.18, text=f"L_ss = {L_ss:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=-0.03*L_ss, y=D/2, text=f"D={D:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False, textangle=-90)
    fig.update_layout(**THEME, title="Vessel Schematic (Side View)",
                      height=260, margin=dict(l=60, r=20, t=40, b=40), showlegend=False)
    fig.update_xaxes(visible=False, range=[-0.08*L_ss, 1.05*L_ss])
    fig.update_yaxes(visible=False, range=[-D*0.35, D*1.15], scaleanchor="x")
    return fig


def plot_schematic_vertical(D, H_total, H_liq, H_gas, n_phases=2,
                             H_water=0, H_oil=0, H_mist=0.15):
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, y0=0, x1=D, y1=H_total,
                  line=dict(color=AMBER, width=2), fillcolor="#21262d")
    if n_phases == 2:
        fig.add_shape(type="rect", x0=0.05*D, y0=0, x1=0.95*D, y1=H_liq,
                      fillcolor="rgba(56,139,253,0.5)", line=dict(width=0))
        fig.add_annotation(x=D/2, y=H_liq/2, text="Liquid",
                           font=dict(color=BLUE), showarrow=False)
        fig.add_annotation(x=D/2, y=H_liq + H_gas/2, text="Gas",
                           font=dict(color=GREY), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0.05*D, y0=0, x1=0.95*D, y1=H_water,
                      fillcolor="rgba(63,185,80,0.5)", line=dict(width=0))
        fig.add_shape(type="rect", x0=0.05*D, y0=H_water, x1=0.95*D, y1=H_water+H_oil,
                      fillcolor="rgba(56,139,253,0.5)", line=dict(width=0))
        for y, text, col in [(H_water/2, "Water", GREEN),
                              (H_water+H_oil/2, "Oil", BLUE) if H_oil > 0.05 else (None, None, None),
                              (H_water+H_oil+H_gas/2, "Gas", GREY)]:
            if y is not None:
                fig.add_annotation(x=D/2, y=y, text=text,
                                   font=dict(color=col), showarrow=False)
    # Demister pad
    H_mist_bot = H_liq + H_gas if n_phases == 2 else H_water + H_oil + H_gas
    if H_mist > 0:
        fig.add_shape(type="rect", x0=0.05*D, y0=H_mist_bot, x1=0.95*D,
                      y1=H_mist_bot + H_mist,
                      fillcolor="rgba(240,165,0,0.35)", line=dict(color=AMBER, width=1))
        fig.add_annotation(x=D*1.05, y=H_mist_bot + H_mist/2, text="Demister",
                           font=dict(color=AMBER, size=8), showarrow=False)
    fig.add_annotation(x=D*1.2, y=H_total/2, text=f"H={H_total:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=D/2, y=-H_total*0.06, text=f"D={D:.2f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.update_layout(**THEME, title="Vessel Schematic (Front View)",
                      height=380, margin=dict(l=40, r=80, t=40, b=40), showlegend=False)
    fig.update_xaxes(visible=False, range=[-0.1*D, 1.5*D])
    fig.update_yaxes(visible=False, range=[-H_total*0.12, H_total*1.1], scaleanchor="x")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# NOZZLE LAYOUT SCHEMATICS (Nozzles & Weir tab)
# ─────────────────────────────────────────────────────────────────────────────

def plot_horizontal_with_nozzles(D, L_ss, nz_pos, weir, is_three):
    fig = go.Figure()
    fig.add_shape(type="rect", x0=0, y0=0, x1=L_ss, y1=D,
                  line=dict(color=AMBER, width=2.5), fillcolor="#1a1f2e")

    if is_three:
        for y0, y1, col in [(0, D*0.25, "rgba(63,185,80,0.25)"),
                             (D*0.25, D*0.60, "rgba(56,139,253,0.20)")]:
            fig.add_shape(type="rect", x0=0, y0=y0, x1=L_ss, y1=y1,
                          line=dict(width=0), fillcolor=col)
        for y, text, col in [(D*0.12, "Water", GREEN), (D*0.42, "Oil", BLUE),
                              (D*0.80, "Gas", GREY)]:
            fig.add_annotation(x=L_ss*0.12, y=y, text=text,
                               font=dict(color=col, size=10), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0, y0=0, x1=L_ss, y1=D*0.50,
                      line=dict(width=0), fillcolor="rgba(56,139,253,0.25)")
        fig.add_annotation(x=L_ss*0.12, y=D*0.25, text="Liquid",
                           font=dict(color=BLUE, size=10), showarrow=False)

    if is_three and weir:
        xw = weir["x_weir_m"]
        for x, h, col, label in [
            (xw,         weir["h_water_weir_m"], GREEN, f"Water weir\n{weir['h_water_weir_m']:.3f} m"),
            (xw+D*0.05, weir["h_oil_weir_m"],   BLUE,  f"Oil weir\n{weir['h_oil_weir_m']:.3f} m"),
        ]:
            fig.add_shape(type="line", x0=x, y0=0, x1=x, y1=h,
                          line=dict(color=col, width=2.5))
            fig.add_annotation(x=x, y=h+D*0.03, text=label,
                               font=dict(color=col, size=9), showarrow=False)

    for name, (px, py, _) in nz_pos.items():
        col = NOZZLE_COLORS.get(name, "#ffffff")
        sym = NOZZLE_SYMBOLS_H.get(name, "circle")
        fig.add_trace(go.Scatter(
            x=[px], y=[py], mode="markers+text",
            marker=dict(symbol=sym, size=12, color=col,
                        line=dict(color="#0d1117", width=1)),
            text=[name], textposition="top center",
            textfont=dict(color=col, size=9), name=name, showlegend=True,
        ))
        if py > D * 0.85:
            fig.add_shape(type="line", x0=px, y0=D, x1=px, y1=D+D*0.12,
                          line=dict(color=col, width=1, dash="dot"))
        elif py < D * 0.15:
            fig.add_shape(type="line", x0=px, y0=0, x1=px, y1=-D*0.12,
                          line=dict(color=col, width=1, dash="dot"))

    fig.add_annotation(x=L_ss/2, y=-D*0.28, text=f"L_ss = {L_ss:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=-L_ss*0.04, y=D/2, text=f"D = {D:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False, textangle=-90)
    fig.update_layout(**THEME, title="Horizontal Separator — Nozzle Layout",
                      height=420, margin=dict(l=70, r=30, t=50, b=60),
                      legend=dict(orientation="h", y=-0.25, font=dict(size=9),
                                  bgcolor="rgba(0,0,0,0)"))
    fig.update_xaxes(visible=False, range=[-L_ss*0.08, L_ss*1.06])
    fig.update_yaxes(visible=False, range=[-D*0.45, D*1.25], scaleanchor="x")
    return fig


def plot_vertical_with_nozzles(D, H_total, nz_pos, best, is_three, me_thickness=0.15):
    fig = go.Figure()
    H_liq = best.get("H_liq_m", best.get("H_oil_m", 0) + best.get("H_water_m", 0))

    fig.add_shape(type="rect", x0=0, y0=0, x1=D, y1=H_total,
                  line=dict(color=AMBER, width=2.5), fillcolor="#1a1f2e")

    if is_three:
        H_w = best.get("H_water_m", 0)
        H_o = best.get("H_oil_m", 0)
        for y0, y1, col in [(0, H_w, "rgba(63,185,80,0.30)"),
                             (H_w, H_w+H_o, "rgba(56,139,253,0.25)")]:
            fig.add_shape(type="rect", x0=0.05*D, y0=y0, x1=0.95*D, y1=y1,
                          fillcolor=col, line=dict(width=0))
        H_g = best.get("H_gas_m", 0.6)
        for y, text, col in [(H_w/2, "Water", GREEN),
                              (H_w+H_o/2, "Oil", BLUE) if H_o > 0.05 else (None,None,None),
                              (H_w+H_o+H_g/2, "Gas", GREY)]:
            if y is not None:
                fig.add_annotation(x=D/2, y=y, text=text,
                                   font=dict(color=col, size=9), showarrow=False)
    else:
        fig.add_shape(type="rect", x0=0.05*D, y0=0, x1=0.95*D, y1=H_liq,
                      fillcolor="rgba(56,139,253,0.25)", line=dict(width=0))
        fig.add_annotation(x=D/2, y=H_liq/2, text="Liquid",
                           font=dict(color=BLUE, size=9), showarrow=False)

    H_mist_bot = H_liq + best.get("H_gas_m", max(0.6, D*0.5))
    H_mist_h   = best.get("H_mist_m", me_thickness)
    if H_mist_h > 0:
        fig.add_shape(type="rect", x0=0.05*D, y0=H_mist_bot,
                      x1=0.95*D, y1=H_mist_bot+H_mist_h,
                      fillcolor="rgba(240,165,0,0.35)", line=dict(color=AMBER, width=1))
        fig.add_annotation(x=D*1.05, y=H_mist_bot+H_mist_h/2, text="Demister",
                           font=dict(color=AMBER, size=8), showarrow=False)

    for name, (px, pz, _) in nz_pos.items():
        col   = NOZZLE_COLORS.get(name, "#ffffff")
        is_top    = pz >= H_total * 0.9
        is_bottom = pz <= 0.01
        y_sym = H_total if is_top else (0 if is_bottom else pz)

        if is_top:
            fig.add_shape(type="line", x0=D/2, y0=H_total, x1=D/2, y1=H_total+D*0.15,
                          line=dict(color=col, width=1.5, dash="dot"))
        elif is_bottom:
            fig.add_shape(type="line", x0=D/2, y0=0, x1=D/2, y1=-D*0.15,
                          line=dict(color=col, width=1.5, dash="dot"))
        else:
            fig.add_shape(type="line", x0=D, y0=pz, x1=D+D*0.2, y1=pz,
                          line=dict(color=col, width=1.5, dash="dot"))

        fig.add_trace(go.Scatter(
            x=[D*0.5], y=[y_sym], mode="markers+text",
            marker=dict(symbol="circle", size=10, color=col,
                        line=dict(color="#0d1117", width=1)),
            text=[name], textposition="middle right",
            textfont=dict(color=col, size=9), name=name, showlegend=True,
        ))

    fig.add_annotation(x=D*1.35, y=H_total/2, text=f"H = {H_total:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.add_annotation(x=D/2, y=-D*0.25, text=f"D = {D:.3f} m",
                       font=dict(color="#e6edf3", size=11), showarrow=False)
    fig.update_layout(**THEME, title="Vertical Separator — Nozzle Layout",
                      height=520, margin=dict(l=40, r=100, t=50, b=60),
                      legend=dict(orientation="v", x=1.05, font=dict(size=9),
                                  bgcolor="rgba(0,0,0,0)"))
    fig.update_xaxes(visible=False, range=[-D*0.3, D*1.8])
    fig.update_yaxes(visible=False, range=[-H_total*0.15, H_total*1.15], scaleanchor="x")
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# MIST EXTRACTOR ΔP CHART
# ─────────────────────────────────────────────────────────────────────────────

def plot_me_dp_curve(me_result):
    """Pressure drop vs gas velocity chart for the selected mist extractor."""
    import numpy as np
    v_max = me_result["v_allow_ms"]
    v_op  = me_result["v_me_ms"]
    dP_c  = me_result["dP_clean_Pa"]
    dP_l  = me_result["dP_loaded_Pa"]

    v_range = np.linspace(0.01, v_max * 1.2, 80)
    dP_range = dP_c * (v_range / v_op)**2

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=v_range, y=dP_range, mode="lines",
        line=dict(color=BLUE, width=2), name="Est. ΔP",
    ))
    fig.add_vline(x=v_op, line=dict(color=AMBER, dash="dash", width=1.5),
                  annotation_text=f"Operating {v_op:.3f} m/s",
                  annotation_font=dict(color=AMBER, size=10))
    fig.add_vline(x=v_max, line=dict(color=RED, dash="dot", width=1.5),
                  annotation_text=f"Max {v_max:.3f} m/s",
                  annotation_font=dict(color=RED, size=10))
    fig.add_hrect(y0=dP_c, y1=dP_l,
                  fillcolor="rgba(240,165,0,0.08)",
                  line=dict(color="rgba(240,165,0,0.3)"),
                  annotation_text="Clean→Loaded range",
                  annotation_font=dict(color=AMBER, size=9))
    fig.update_layout(**THEME, **GRID,
                      title="Mist Extractor ΔP vs Gas Velocity (estimated)",
                      xaxis_title="Gas Velocity [m/s]",
                      yaxis_title="Pressure Drop [Pa]",
                      height=300, margin=dict(l=60, r=20, t=45, b=40))
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# SIZING TABLE BAR CHART
# ─────────────────────────────────────────────────────────────────────────────

def plot_length_vs_diameter(df, best):
    bar_y = "L_eff_m" if "L_eff_m" in df.columns else "H_total_m"
    label = "Length [m]" if "L_eff_m" in df.columns else "Height [m]"
    fig = px.bar(df, x="D_m", y=bar_y,
                 color="D_m",
                 color_continuous_scale=[[0, "#21262d"], [1, AMBER]],
                 labels={"D_m": "Diameter [m]", bar_y: label},
                 title="Effective Length / Height vs Diameter")
    fig.add_vline(x=best["D_m"], line=dict(color=AMBER, dash="dash"),
                  annotation_text="Selected")
    fig.update_layout(**THEME, height=300, margin=dict(l=50, r=20, t=40, b=40))
    return fig
