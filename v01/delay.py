#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Auswertung und Plot der Verzoegerungs-Koinzidenzmessung (delay.csv).

Rekonstruiert aus den Spalten "Count left" / "Count right" die symmetrische
T_VZ-Kurve, bestimmt das Plateau (Mittelwert des flachen Bereichs) und die
Halbwertsbreite (FWHM) ueber die Halb-Maximum-Schnittpunkte und stellt das
Ganze im Stil der Vorlage dar.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# EINSTELLUNGEN (bei Bedarf anpassen)
# ----------------------------------------------------------------------
CSV_PFAD       = "raw/delay.csv"
PLATEAU_FENSTER = 12      # Punkte mit |T_VZ| <= diesem Wert [ns] bilden das Plateau
TAU_REF        = 20.0     # Referenzzeit fuer Delta t_K = 2*TAU_REF - FWHM  [ns]
MESSZEIT_S     = 20       # nur zur Info (Counts "in 20s")
OUTPUT_PNG     = "delay_plot.png"

# ----------------------------------------------------------------------
# DATEN EINLESEN UND T_VZ-ACHSE REKONSTRUIEREN
# ----------------------------------------------------------------------
df = pd.read_csv("raw/delay.csv")
dt    = df["Delta t[ns]"].to_numpy(dtype=float)
left  = df["Count left (in 20s)"].to_numpy(dtype=float)
right = df["Count right (in 20s)"].to_numpy(dtype=float)

# left -> negative Achse, right -> positive Achse, Delta t = 0 -> T_VZ = 0
t_neg = -dt[dt > 0]
c_neg = left[dt > 0]

t_pos = dt[dt > 0]
c_pos = right[dt > 0]

# Nullpunkt (Delta t = 0): nutzt den vorhandenen left-Wert
t_zero = np.array([0.0])
c_zero = np.array([left[dt == 0][0]])

# entferne fehlende (NaN) Eintraege
mask_neg = ~np.isnan(c_neg)
mask_pos = ~np.isnan(c_pos)

t = np.concatenate([t_neg[mask_neg], t_zero, t_pos[mask_pos]])
c = np.concatenate([c_neg[mask_neg], c_zero, c_pos[mask_pos]])

# nach T_VZ sortieren
order = np.argsort(t)
t = t[order]
c = c[order]

# Poisson-Fehler
c_err = np.sqrt(c)

# ----------------------------------------------------------------------
# PLATEAU BESTIMMEN (Mittelwert des flachen Zentralbereichs)
# ----------------------------------------------------------------------
plateau_mask = np.abs(t) <= PLATEAU_FENSTER
plateau      = c[plateau_mask].mean()
plateau_err  = c[plateau_mask].std(ddof=1) / np.sqrt(plateau_mask.sum())
half_max     = plateau / 2.0

# ----------------------------------------------------------------------
# HALBWERTSBREITE: SCHNITTPUNKTE MIT half_max (lineare Interpolation)
# ----------------------------------------------------------------------
def schnittpunkt(t_arr, c_arr, level, seite):
    """Findet von aussen nach innen / innen nach aussen den Durchgang durch 'level'."""
    if seite == "links":
        idx = np.where(t_arr < 0)[0]
    else:
        idx = np.where(t_arr > 0)[0]
    ti, ci = t_arr[idx], c_arr[idx]
    o = np.argsort(ti)
    ti, ci = ti[o], ci[o]
    # benachbarte Paare suchen, bei denen level ueberschritten wird
    for k in range(len(ti) - 1):
        a, b = ci[k], ci[k + 1]
        if (a - level) * (b - level) <= 0 and a != b:
            frac = (level - a) / (b - a)
            return ti[k] + frac * (ti[k + 1] - ti[k])
    return np.nan

t_links  = schnittpunkt(t, c, half_max, "links")
t_rechts = schnittpunkt(t, c, half_max, "rechts")
fwhm     = t_rechts - t_links
dt_K     = 2 * TAU_REF - fwhm

print(f"Plateau         = {plateau:.1f} +/- {plateau_err:.1f} counts")
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

# Messwerte mit Fehlerbalken
ax.errorbar(t, c, yerr=c_err, fmt="+", color="black", ecolor="black",
            elinewidth=0.9, capsize=2.5, markersize=8, label="Messwerte")

# Plateau-Linie (rot, durchgezogen) ueber den Plateaubereich
t_plat = t[plateau_mask]
ax.plot([t_plat.min(), t_plat.max()], [plateau, plateau],
        color="red", lw=1.3, label="Plateau")

# Halbwertsbreite (rot, gestrichelt) zwischen den Schnittpunkten
ax.plot([t_links, t_rechts], [half_max, half_max],
        color="red", lw=1.2, ls="--", label="Halbwertsbreite")

# senkrechte gestrichelte Linien an den Schnittpunkten
for tc in (t_links, t_rechts):
    ax.axvline(tc, color="black", lw=0.9, ls="--")
    ax.annotate(f"{tc:.0f}", xy=(tc, ax.get_ylim()[1]),
                xytext=(tc, plateau * 1.45), color="red",
                ha="center", fontsize=12, fontweight="bold")

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