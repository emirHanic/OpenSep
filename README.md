![image](./LOGO.png)

# Separator Sizing Calculator — API 12J / GPSA EDB

A Streamlit web application for preliminary sizing of oil & gas phase separators,
strictly following **API 12J** (Specification for Oil and Gas Separators) and
the **GPSA Engineering Data Book**, in SI units throughout.

## Capabilities

| Configuration                          | Supported |
| -------------------------------------- | --------- |
| Horizontal Two-Phase (Gas/Liquid)      | ✅        |
| Horizontal Three-Phase (Gas/Oil/Water) | ✅        |
| Vertical Two-Phase (Gas/Liquid)        | ✅        |
| Vertical Three-Phase (Gas/Oil/Water)   | ✅        |

## Methods & Standards

- **Droplet settling**: Iterative drag coefficient (Schiller-Naumann) per GPSA Ch. 7
- **Max gas velocity**: Souders-Brown equation with pressure-dependent K factor (API 12J Table 1 / GPSA Fig. 7-3)
- **Liquid retention**: API 12J Table 2 — based on oil API gravity and operating temperature
- **L/D selection**: API 12J preferred range 3–5 (horizontal), 1.5–4 (vertical)
- **Gas density**: Real gas law with Z-factor correction
- **Three-phase separation**: Oil/water interface sizing using droplet rise/settling velocities

## Installation

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Notes

- All inputs and outputs are in SI units (m, kg, s, K, kPa)
- Results are for **preliminary sizing only**
- Final design must comply with API 12J, ASME VIII, and applicable local codes
- Review by a qualified process/mechanical engineer is required
