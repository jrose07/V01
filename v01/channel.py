#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zeitkalibrierung des MCA (channel.csv).

Bekannte Zeitdifferenzen t [us] werden gegen die gemessene Kanalnummer
aufgetragen. Es wird eine lineare Ausgleichsrechnung
    Delta_t(K) = a * K + b
mit scipy.optimize.curve_fit durchgefuehrt und im Stil der Vorlage geplottet.
"""
import uncertainties.unumpy as unp
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from addons import write, add, latex_float, tab_to_latex as tab2tex
from scipy.stats import linregress
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
import scipy.constants as const
import matplotlib.pyplot as plt

dir     = "content/plots/"
dir_tab = "content/tables/"

# ----------------------------------------------------------------------
# EINSTELLUNGEN
# ----------------------------------------------------------------------
CSV_PFAD   = "raw/channel.csv"
OUTPUT_PNG = dir + "channel_plot.png"


# ----------------------------------------------------------------------
# HILFSFUNKTION: deutsche Zahlen robust einlesen
# ----------------------------------------------------------------------
def zu_float(x):
    """Wandelt Strings mit Dezimalkomma in float um; korrigiert fehlende Kommas."""
    if pd.isna(x):
        return np.nan
    s = str(x).strip().replace(",", ".")
    if s == "":
        return np.nan
    v = float(s)
    # Artefakt: Komma fehlte beim Export (z.B. "246026" -> 246.026)
    if v > 1000:
        v = v / 1000.0
    return v


# ----------------------------------------------------------------------
# DATEN EINLESEN
# ----------------------------------------------------------------------
df = pd.read_csv(CSV_PFAD)
df.columns    = ["t_us", "channel"]
df["t_us"]    = df["t_us"].map(zu_float)
df["channel"] = df["channel"].map(zu_float)
df = df.dropna().reset_index(drop=True)

channel = df["channel"].to_numpy()
t_us    = df["t_us"].to_numpy()


# ----------------------------------------------------------------------
# LINEARE AUSGLEICHSRECHNUNG  Delta_t = a*K + b
# ----------------------------------------------------------------------
def gerade(K, a, b):
    return a * K + b

popt, pcov = curve_fit(gerade, channel, t_us)
a, b   = popt
da, db = np.sqrt(np.diag(pcov))

# als ufloat fuer schoenere Ausgabe / Weiterverwendung
a_u = ufloat(a, da)
b_u = ufloat(b, db)

# Bestimmtheitsmass R^2
resid = t_us - gerade(channel, *popt)
r2 = 1 - np.sum(resid**2) / np.sum((t_us - t_us.mean())**2)

print(f"Steigung          a = {a_u*1000} ns/Kanal")
print(f"                    = {a_u} us/Kanal")
print(f"Achsenabschnitt   b = {b_u} us")
print(f"R^2                 = {r2:.6f}")


# ----------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------
plt.rcParams.update({"font.size": 12, "font.family": "serif"})
fig, ax = plt.subplots(figsize=(8, 6))

# Ausgleichsgerade (etwas ueber den Datenbereich hinaus)
K_fit = np.linspace(channel.min() - 20, channel.max() + 20, 200)
ax.plot(K_fit, gerade(K_fit, *popt), color="black", lw=1.2,
        label="Ausgleichsfunktion", zorder=1)

# Messwerte
ax.plot(channel, t_us, "x", color="red", markersize=7, mew=1.4,
        label="Messwerte", zorder=2)

ax.set_xlabel("Kanalnummer")
ax.set_ylabel(r"$\Delta t$ / µs")
ax.grid(True, color="0.85", lw=0.6)
ax.legend(loc="upper left", framealpha=1.0, edgecolor="black")

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=150)
print(f"\nGespeichert: {OUTPUT_PNG}")