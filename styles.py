"""
styles.py — OpenSep
All CSS and HTML helper strings for the Streamlit UI.
Centralised here so visual changes never touch business logic.
"""

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@300;400;600;700&display=swap');

html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; }

.stApp { background-color: #0d1117; color: #e6edf3; }

[data-testid="stSidebar"] {
    background-color: #161b22;
    border-right: 1px solid #30363d;
}
[data-testid="stSidebar"] .stMarkdown h2,
[data-testid="stSidebar"] .stMarkdown h3 {
    color: #f0a500; font-family: 'IBM Plex Mono', monospace;
    font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.15em;
    border-bottom: 1px solid #30363d; padding-bottom: 0.3rem; margin-top: 1.2rem;
}

.result-card {
    background: #161b22; border: 1px solid #30363d;
    border-left: 3px solid #f0a500; border-radius: 6px;
    padding: 1rem 1.25rem; margin: 0.4rem 0;
}
.result-card .label {
    font-size: 0.72rem; color: #8b949e; text-transform: uppercase;
    letter-spacing: 0.1em; font-family: 'IBM Plex Mono', monospace;
}
.result-card .value {
    font-size: 1.5rem; font-weight: 700; color: #f0a500;
    font-family: 'IBM Plex Mono', monospace;
}
.result-card .unit { font-size: 0.85rem; color: #8b949e; margin-left: 0.3rem; }

.warn-card {
    background: #1c1810; border: 1px solid #f0a500; border-radius: 6px;
    padding: 0.75rem 1rem; margin: 0.4rem 0; color: #f0a500; font-size: 0.85rem;
}
.ok-card {
    background: #0d1f12; border: 1px solid #238636; border-radius: 6px;
    padding: 0.75rem 1rem; margin: 0.4rem 0; color: #3fb950; font-size: 0.85rem;
}

.stTabs [data-baseweb="tab-list"] {
    background-color: #161b22; border-bottom: 1px solid #30363d; gap: 0;
}
.stTabs [data-baseweb="tab"] {
    background-color: transparent; color: #8b949e;
    border-bottom: 2px solid transparent; padding: 0.6rem 1.2rem;
    font-family: 'IBM Plex Mono', monospace; font-size: 0.8rem;
    text-transform: uppercase; letter-spacing: 0.08em;
}
.stTabs [aria-selected="true"] {
    color: #f0a500 !important; border-bottom: 2px solid #f0a500 !important;
    background-color: transparent !important;
}

h1 { font-family: 'IBM Plex Mono', monospace !important; color: #f0a500 !important;
     font-size: 1.6rem !important; letter-spacing: -0.02em; }
h2, h3 { font-family: 'IBM Plex Mono', monospace !important; color: #e6edf3 !important; }

.stNumberInput input, .stSelectbox select, .stSlider {
    background-color: #21262d !important; border-color: #30363d !important;
    color: #e6edf3 !important;
}

[data-testid="metric-container"] {
    background: #161b22; border: 1px solid #30363d; border-radius: 6px; padding: 0.8rem 1rem;
}
[data-testid="metric-container"] label {
    color: #8b949e !important; font-size: 0.72rem !important;
    font-family: 'IBM Plex Mono', monospace !important;
    text-transform: uppercase; letter-spacing: 0.1em;
}
[data-testid="metric-container"] [data-testid="metric-value"] {
    color: #f0a500 !important; font-family: 'IBM Plex Mono', monospace !important;
}

hr { border-color: #30363d !important; }
.stInfo { background-color: #1c2128; border-color: #388bfd; }

.section-pill {
    display: inline-block; background: #f0a50020; border: 1px solid #f0a50060;
    color: #f0a500; font-family: 'IBM Plex Mono', monospace; font-size: 0.68rem;
    text-transform: uppercase; letter-spacing: 0.15em; padding: 0.15rem 0.6rem;
    border-radius: 100px; margin-bottom: 0.6rem;
}
</style>
"""


def result_card(label, value, unit=""):
    return (
        f"<div class='result-card'>"
        f"<div class='label'>{label}</div>"
        f"<div class='value'>{value}<span class='unit'>{unit}</span></div>"
        f"</div>"
    )


def warn_card(text):
    return f"<div class='warn-card'>{text}</div>"


def ok_card(text):
    return f"<div class='ok-card'>{text}</div>"


def section_pill(text):
    return f"<div class='section-pill'>{text}</div>"


HEADER_HTML = """
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
"""

FOOTER_HTML = """
<hr>
<p style='text-align:center;color:#8b949e;font-size:0.72rem;font-family:IBM Plex Mono'>
API 12J · GPSA Engineering Data Book · SI Units Throughout · For Preliminary Sizing Only
</p>
"""

THEORY_MD = """
## Theory & Reference Equations

### 1. Fluid Properties

**Gas Density** (Real Gas Law):
$$\\rho_g = \\frac{P \\cdot M_w}{Z \\cdot R \\cdot T}$$

---

### 2. Droplet Settling Velocity

**Stokes' Law** (Re < 0.4):
$$v_t = \\frac{d_p^2 (\\rho_L - \\rho_G) g}{18 \\mu_G}$$

**Schiller-Naumann drag** (API 12J §5 / GPSA Ch. 7):
$$v_t = \\sqrt{\\frac{4 d_p (\\rho_L - \\rho_G) g}{3 C_D \\rho_G}}$$

---

### 3. Souders-Brown (API 12J Table 1 / GPSA Fig. 7-3)

$$v_{SB} = K \\sqrt{\\frac{\\rho_L - \\rho_G}{\\rho_G}}$$

---

### 4. Horizontal Two-Phase (API 12J §5.3)

$$L_{eff} \\geq \\frac{v_G \\cdot D}{v_t}, \\quad L_{eff} \\geq \\frac{Q_L \\cdot t_R}{A_L}$$
$$L_{ss} = L_{eff} + 0.3D$$

---

### 5. Vertical Two-Phase (API 12J §5.2)

$$H_T = H_{liq} + H_{gas} + H_{mist} + H_{inlet}$$
$H_{gas} \\geq 0.6$ m with ME, $\\geq 0.9$ m without.

---

### 6. Three-Phase Sizing (API 12J §5.4)

Oil droplet rise in water: $v_{t,o/w} = d_{p,o}^2 (\\rho_W - \\rho_O) g / (18 \\mu_W)$

Water droplet settling in oil: $v_{t,w/o} = d_{p,w}^2 (\\rho_W - \\rho_O) g / (18 \\mu_O)$

---

### 7. Design Criteria (API 12J / GPSA)

| Parameter | Horizontal | Vertical |
|-----------|-----------|---------|
| Preferred L/D or H/D | 3 – 5 | 1.5 – 4 |
| Max gas velocity | ≤ $v_{SB}$ | ≤ $\\min(v_{SB}, 0.75v_t)$ |

---

### 8. Mist Extractor (API 12J §8)

$$A_{ME} = \\frac{Q_G}{K_{ME}\\sqrt{(\\rho_L-\\rho_G)/\\rho_G}}, \\quad
\\Delta P = \\Delta P_{clean}\\left(\\frac{v_{ME}}{v_{ref}}\\right)^2$$

---

### 9. Inlet Device (API 12J §6.3)

$K_{eff} = K_{base} \\times f_K$, $\\quad L_{eff} \\times f_L$, $\\quad Q_{eff} = Q(1-f_{pre-sep})$

---

### 10. Nozzle Sizing (API 12J §6)

Inlet: $\\rho v^2 \\leq 1490\\ \\text{kg/(m·s}^2)$.  
Gas outlet: $v \\leq 18$ m/s.  Liquid outlet: $v \\leq 3$ m/s.  
Next standard DN per ASME B36.10.

---

### 11. Weir Design (API 12J §5.4.3)

$h_{water\\ weir} = 0.25D$, $\\quad h_{oil\\ weir} = 0.60D$,
$\\quad \\Delta h_{oil} = h_{oil}-h_{water}$

---

### References
- **API 12J** — Specification for Oil and Gas Separators, 8th Edition
- **GPSA Engineering Data Book** — Chapter 7
- Schiller & Naumann (1935); Souders & Brown (1934)
"""
