#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Langzeitmessung zur Lebensdauer der Myonen (.Spe-Datei des MCA).

Liest das Maestro/Ortec-.Spe-File aus (Messzeit, Startzeit -> Endzeit,
registrierte Stoppimpulse), rechnet die Kanalnummer mit der Kalibrierung
aus channel.py in eine Zeit um und fuehrt eine exponentielle
Ausgleichsrechnung
    N(t) = N_0 * exp(-lambda * t) + U
mit scipy.optimize.curve_fit durch. Plot im Stil der Vorlage (Abb. 4).
"""
import uncertainties.unumpy as unp
import numpy as np
import datetime as dt
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
SPE_PFAD   = "raw/2026_0611.Spe"
OUTPUT_PNG = dir + "lifetime_plot.png"

# Zeitkalibrierung aus channel.py:  Delta_t = CAL_A * Kanal + CAL_B
CAL_A = 0.021685      # us / Kanal
CAL_B = 0.155200      # us
SKIP_FIRST = 3        # erste Kanaele nicht in den Fit einbeziehen


# ----------------------------------------------------------------------
# .Spe-DATEI EINLESEN
# ----------------------------------------------------------------------
def lies_spe(pfad):
    with open(pfad) as f:
        lines = [l.rstrip("\r\n") for l in f]

    def sektion(name):
        i = lines.index(name)
        out = []
        for l in lines[i + 1:]:
            if l.startswith("$"):
                break
            out.append(l)
        return out

    start    = dt.datetime.strptime(sektion("$DATE_MEA:")[0], "%m/%d/%Y %H:%M:%S")
    live, real = (int(x) for x in sektion("$MEAS_TIM:")[0].split())
    ch0, ch1   = (int(x) for x in sektion("$DATA:")[0].split())
    n_ch       = ch1 - ch0 + 1
    di         = lines.index("$DATA:")
    counts     = np.array([int(x) for x in lines[di + 2: di + 2 + n_ch]], dtype=float)
    return start, live, real, counts


start, live, real, counts = lies_spe(SPE_PFAD)
ende = start + dt.timedelta(seconds=real)

kanal = np.arange(len(counts))
t_us  = CAL_A * kanal + CAL_B

N_stop_spektrum = int(counts.sum())          # registrierte Stoppimpulse im Spektrum
letzter_kanal   = int(np.max(np.where(counts > 0)))

print("----- Metadaten aus der .Spe-Datei -----")
print(f"Startzeit        : {start:%d.%m.%Y %H:%M:%S}")
print(f"Endzeit          : {ende:%d.%m.%Y %H:%M:%S}")
print(f"Messdauer (real) : {real} s = {real/3600:.2f} h = {real/86400:.2f} Tage")
print(f"Messdauer (live) : {live} s")
print(f"N_stop (Spektrum): {N_stop_spektrum}")
print(f"letzter Kanal >0 : {letzter_kanal}")
print("(N_start ist nicht in der Datei gespeichert -> Zaehleranzeige verwenden)\n")


# ----------------------------------------------------------------------
# MASKEN: einbezogene / nicht einbezogene Messwerte
# ----------------------------------------------------------------------
einbezogen = (kanal >= SKIP_FIRST) & (kanal <= letzter_kanal)
ausgelassen = ~einbezogen


# ----------------------------------------------------------------------
# AUSGLEICHSRECHNUNG  N(t) = N0 * exp(-lambda t) + U
# ----------------------------------------------------------------------
def zerfall(t, N0, lam, U):
    return N0 * np.exp(-lam * t) + U

p0 = [counts.max(), 0.46, 0.3]
popt, pcov = curve_fit(zerfall, t_us[einbezogen], counts[einbezogen],
                       p0=p0, maxfev=10000)
perr = np.sqrt(np.diag(pcov))

N0_u  = ufloat(popt[0], perr[0])
lam_u = ufloat(popt[1], perr[1])
U_u   = ufloat(popt[2], perr[2])
tau_u = 1 / lam_u

print("----- Ausgleichsrechnung -----")
print(f"N_0    = {N0_u}")
print(f"lambda = {lam_u} 1/us")
print(f"U      = {U_u}")
print(f"tau    = {tau_u} us")


# ----------------------------------------------------------------------
# PLOT (Stil Abbildung 4)
# ----------------------------------------------------------------------
plt.rcParams.update({"font.size": 12, "font.family": "serif"})
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(t_us[einbezogen], counts[einbezogen], ".", color="red", ms=4,
        label="Messwerte")
ax.plot(t_us[ausgelassen], counts[ausgelassen], ".", color="blue", ms=4,
        label="nicht einbezogene Messwerte")

t_fit = np.linspace(t_us[einbezogen].min(), t_us[einbezogen].max(), 500)
ax.plot(t_fit, zerfall(t_fit, *popt), color="black", lw=1.0,
        label="Ausgleichsfunktion")

ax.set_xlabel(r"$t$ / µs")
ax.set_ylabel("Counts")
ax.grid(True, color="0.9", lw=0.6)
ax.legend(loc="upper right", framealpha=1.0, edgecolor="black")

fig.tight_layout()
fig.savefig(OUTPUT_PNG, dpi=150)
print(f"\nGespeichert: {OUTPUT_PNG}")