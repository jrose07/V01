#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auswertung und Plot der Verzoegerungs-Koinzidenzmessung (delay.csv).

Rekonstruiert aus den Spalten "Count left" / "Count right" die symmetrische
T_VZ-Kurve und bestimmt das Plateau (Mittelwert des flachen Bereichs).
Die Halbwertsbreite wird ueber eine Gerade durch die ZWEI das Halb-Maximum
einrahmenden Punkte bestimmt (einer knapp darunter, einer knapp darueber).
Der Schnittpunkt dieser Geraden mit der Halb-Maximum-Linie liefert die linke
bzw. rechte Flankenzeit; die Unsicherheit folgt aus den Poisson-Fehlern der
beiden Counts.
"""
import os
import uncertainties.unumpy as unp
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from addons import write, add, latex_float, tab_to_latex as tab2tex
from scipy.stats import linregress
from uncertainties import ufloat, correlated_values
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
import scipy.constants as const
import matplotlib.pyplot as plt

dir     = "content/plots/"
dir_tab = "content/tables/"
os.makedirs(dir, exist_ok=True)

# ----------------------------------------------------------------------
# EINSTELLUNGEN
# ----------------------------------------------------------------------
CSV_PFAD        = "raw/delay.csv"
OUTPUT_PNG      = dir + "delay_plot.png"
PLATEAU_FENSTER = 12      # |T_VZ| <= dieser Wert [ns] bildet das Plateau
TAU_REF         = 20.0    # Pulsbreite; Aufloesung = |2*TAU_REF - FWHM|  [ns]


# ----------------------------------------------------------------------
# DATEN EINLESEN UND T_VZ-ACHSE REKONSTRUIEREN
# ----------------------------------------------------------------------
df = pd.read_csv(CSV_PFAD)
dt    = df["Delta t[ns]"].to_numpy(dtype=float)
left  = df["Count left (in 20s)"].to_numpy(dtype=float)
right = df["Count right (in 20s)"].to_numpy(dtype=float)

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
# PLATEAU (mit Unsicherheit)  ->  Halb-Maximum
# ----------------------------------------------------------------------
plateau_mask = np.abs(t) <= PLATEAU_FENSTER
plateau      = c[plateau_mask].mean()
plateau_err  = c[plateau_mask].std(ddof=1) / np.sqrt(plateau_mask.sum())
plateau_u    = ufloat(plateau, plateau_err)
half_max_u   = plateau_u / 2
half_max     = half_max_u.n


# ----------------------------------------------------------------------
# FLANKEN: Gerade durch die zwei das Halb-Maximum einrahmenden Punkte
# ----------------------------------------------------------------------
def flanken_fit(seite):
    """Gerade durch die zwei benachbarten Punkte, die half_max einrahmen.
    Rueckgabe: (t_cross [ufloat], coeffs [m, b], t_used)."""
    msk = (t < 0) if seite == "links" else (t > 0)
    ti, ci = t[msk], c[msk]
    o = np.argsort(ti)
    ti, ci = ti[o], ci[o]

    # erstes benachbartes Paar, das das Halb-Maximum kreuzt
    k = None
    for j in range(len(ti) - 1):
        if (ci[j] - half_max) * (ci[j + 1] - half_max) <= 0 and ci[j] != ci[j + 1]:
            k = j
            break
    if k is None:
        raise ValueError(f"Keine Flanke gefunden ({seite}): Halb-Maximum nicht gekreuzt.")

    # die zwei Punkte mit Poisson-Fehler sqrt(N)
    t_lo, t_hi = ti[k], ti[k + 1]
    c_lo = ufloat(ci[k],     np.sqrt(ci[k]))
    c_hi = ufloat(ci[k + 1], np.sqrt(ci[k + 1]))

    # Gerade durch zwei Punkte -> Schnittpunkt mit half_max
    m_u = (c_hi - c_lo) / (t_hi - t_lo)
    b_u = c_lo - m_u * t_lo
    t_cross = (half_max_u - b_u) / m_u

    coeffs = [m_u.n, b_u.n]                  # nur fuers Plotten
    return t_cross, coeffs, np.array([t_lo, t_hi])


t_links_u,  coef_l, ti_l = flanken_fit("links")
t_rechts_u, coef_r, ti_r = flanken_fit("rechts")

fwhm_u = t_rechts_u - t_links_u
dt_K_u = 2 * TAU_REF - fwhm_u

print(f"Plateau         = {plateau_u} counts")
print(f"Halb-Maximum    = {half_max:.1f} counts")
print(f"Punkte links    = {ti_l} ns")
print(f"Punkte rechts   = {ti_r} ns")
print(f"Schnittpunkt L  = {t_links_u} ns")
print(f"Schnittpunkt R  = {t_rechts_u} ns")
print(f"Halbwertsbreite = {fwhm_u} ns")
print(f"Aufloesungszeit = |2*{TAU_REF:.0f} - FWHM| = {abs(dt_K_u.n):.2f} +/- {dt_K_u.s:.2f} ns")


# ----------------------------------------------------------------------
# PLOT
# ----------------------------------------------------------------------
plt.rcParams.update({"font.size": 11, "font.family": "serif",
                     "mathtext.fontset": "cm"})
fig, ax = plt.subplots(figsize=(8, 5.2))

# Farbpalette
C_DATA  = "#2B3A55"   # Messwerte
C_PLAT  = "#2A9D8F"   # Plateau
C_FLANK = "#E76F51"   # Flankengeraden
C_HWB   = "#6D4C41"   # Halbwertsbreite

ax.grid(True, color="0.9", lw=0.6, zorder=0)
ax.set_axisbelow(True)

# Messwerte
ax.errorbar(t, c, yerr=c_err, fmt="o", color=C_DATA, ecolor=C_DATA,
            elinewidth=0.8, capsize=2, markersize=3.5, alpha=0.9,
            label="Messwerte", zorder=3)

# Plateau als Linie + Fehlerband
t_plat = t[plateau_mask]
ax.hlines(plateau, t_plat.min(), t_plat.max(), color=C_PLAT, lw=2,
          zorder=4, label="Plateau")

# Halb-Maximum-Niveau + FWHM-Doppelpfeil
ax.hlines(half_max, t_links_u.n, t_rechts_u.n, color=C_HWB, lw=1,
          ls=(0, (4, 3)), zorder=4)
ax.annotate("", xy=(t_rechts_u.n, half_max), xytext=(t_links_u.n, half_max),
            arrowprops=dict(arrowstyle="<->", color=C_HWB, lw=1.3))
ax.text((t_links_u.n + t_rechts_u.n) / 2, half_max + 0.045 * c.max(),
        rf"FWHM $= {fwhm_u.n:.1f}$ ns", color=C_HWB, ha="center",
        fontsize=10, style="italic")

# Flankengeraden: nur zwischen den zwei verwendeten Punkten
for k, (ti_, coef_) in enumerate([(ti_l, coef_l), (ti_r, coef_r)]):
    xx = np.array([ti_[0], ti_[1]])
    ax.plot(xx, coef_[0] * xx + coef_[1], color=C_FLANK, lw=1.8,
            marker="o", markersize=6, mfc="white", mec=C_FLANK, mew=1.5,
            label="Flankengerade" if k == 0 else None, zorder=5)

# Schnittpunkte: dezente Verticals + Wertelabels
for tc in (t_links_u.n, t_rechts_u.n):
    ax.axvline(tc, color="0.6", lw=0.8, ls=":", zorder=1)
    ax.text(tc, 0.045 * c.max(), f"{tc:.1f}", color=C_FLANK, ha="center",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))

ax.set_xlabel(r"$T_{\mathrm{VZ}} \;/\; \mathrm{ns}$")
ax.set_ylabel("Counts")
ax.set_xlim(t.min() - 5, t.max() + 5)
ax.set_ylim(0, c.max() * 1.35)

handles, labels = ax.get_legend_handles_labels()
reihenfolge = ["Messwerte", "Plateau", "Flankengerade"]
idx = [labels.index(l) for l in reihenfolge if l in labels]
ax.legend([handles[i] for i in idx], [labels[i] for i in idx],
          loc="upper right", framealpha=0.95, edgecolor="0.8", fontsize=10)

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=150)
print(f"\nGespeichert: {OUTPUT_PNG}")