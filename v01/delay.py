#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auswertung und Plot der Verzoegerungs-Koinzidenzmessung (delay.csv).

Rekonstruiert aus den Spalten "Count left" / "Count right" die symmetrische
T_VZ-Kurve, bestimmt das Plateau (Mittelwert des flachen Bereichs) und die
Halbwertsbreite (FWHM) ueber die Halb-Maximum-Schnittpunkte und stellt das
Ganze im Stil der Vorlage dar.
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
CSV_PFAD        = "raw/delay.csv"
OUTPUT_PNG      = dir + "delay_plot.png"
PLATEAU_FENSTER = 12      # |T_VZ| <= dieser Wert [ns] bildet das Plateau
TAU_REF         = 20.0    # Referenzzeit fuer Delta t_K = 2*TAU_REF - FWHM  [ns]


# ----------------------------------------------------------------------
# DATEN EINLESEN UND T_VZ-ACHSE REKONSTRUIEREN
# ----------------------------------------------------------------------
df = pd.read_csv(CSV_PFAD)
dt    = df["Delta t[ns]"].to_numpy(dtype=float)
left  = df["Count left (in 20s)"].to_numpy(dtype=float)
right = df["Count right (in 20s)"].to_numpy(dtype=float)

# left -> negative Achse, right -> positive Achse, Delta t = 0 -> T_VZ = 0
t_neg, c_neg = -dt[dt > 0], left[dt > 0]
t_pos, c_pos =  dt[dt > 0], right[dt > 0]
t_zero       = np.array([0.0])
c_zero       = np.array([left[dt == 0][0]])

mask_neg = ~np.isnan(c_neg)
mask_pos = ~np.isnan(c_pos)

t = np.concatenate([t_neg[mask_neg], t_zero, t_pos[mask_pos]])
c = np.concatenate([c_neg[mask_neg], c_zero, c_pos[mask_pos]])

order = np.argsort(t)
t, c = t[order], c[order]
c_err = np.sqrt(c)


# ----------------------------------------------------------------------
# PLATEAU UND HALBWERTSBREITE
# ----------------------------------------------------------------------
plateau_mask = np.abs(t) <= PLATEAU_FENSTER
plateau      = c[plateau_mask].mean()
plateau_err  = c[plateau_mask].std(ddof=1) / np.sqrt(plateau_mask.sum())
plateau_u    = ufloat(plateau, plateau_err)
half_max     = plateau / 2.0


def schnittpunkt(t_arr, c_arr, level, seite):
    """Lineare Interpolation des Durchgangs durch 'level'."""
    idx = np.where(t_arr < 0)[0] if seite == "links" else np.where(t_arr > 0)[0]
    ti, ci = t_arr[idx], c_arr[idx]
    o = np.argsort(ti)
    ti, ci = ti[o], ci[o]
    for k in range(len(ti) - 1):
        a_, b_ = ci[k], ci[k + 1]
        if (a_ - level) * (b_ - level) <= 0 and a_ != b_:
            frac = (level - a_) / (b_ - a_)
            return ti[k] + frac * (ti[k + 1] - ti[k])
    return np.nan


t_links  = schnittpunkt(t, c, half_max, "links")
t_rechts = schnittpunkt(t, c, half_max, "rechts")
fwhm     = t_rechts - t_links
dt_K     = 2 * TAU_REF - fwhm

print(f"Plateau         = {plateau_u} counts")
print(f"Halb-Maximum    = {half_max:.1f} counts")
print(f"Schnittpunkt L  = {t_links:.2f} ns")
print(f"Schnittpunkt R  = {t_rechts:.2f} ns")
print(f"Halbwertsbreite = {fwhm:.2f} ns")
print(f"Delta t_K       = 2*{TAU_REF:.0f} ns - {fwhm:.1f} ns = {dt_K:.2f} ns")


# ----------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------
plt.rcParams.update({"font.size": 12, "font.family": "serif"})
fig, ax = plt.subplots(figsize=(8, 6))

ax.errorbar(t, c, yerr=c_err, fmt="+", color="black", ecolor="black",
            elinewidth=0.9, capsize=2.5, markersize=8, label="Messwerte")

t_plat = t[plateau_mask]
ax.plot([t_plat.min(), t_plat.max()], [plateau, plateau],
        color="red", lw=1.3, label="Plateau")

ax.plot([t_links, t_rechts], [half_max, half_max],
        color="red", lw=1.2, ls="--", label="Halbwertsbreite")

for tc in (t_links, t_rechts):
    ax.axvline(tc, color="black", lw=0.9, ls="--")
    ax.annotate(f"{tc:.0f}", xy=(tc, plateau * 1.45),
                color="red", ha="center", fontsize=12, fontweight="bold")

ax.set_xlabel(r"$T_{\mathrm{VZ}}$ [ns]")
ax.set_ylabel("Counts")
ax.set_xlim(t.min() - 5, t.max() + 5)
ax.set_ylim(0, c.max() * 1.45)

handles, labels = ax.get_legend_handles_labels()
reihenfolge = ["Messwerte", "Plateau", "Halbwertsbreite"]
idx = [labels.index(l) for l in reihenfolge if l in labels]
ax.legend([handles[i] for i in idx], [labels[i] for i in idx],
          loc="upper center", framealpha=1.0, edgecolor="black")

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=150)
print(f"\nGespeichert: {OUTPUT_PNG}")