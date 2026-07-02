#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bestimmung der Untergrundrate (analog Kapitel 3.4 des Altprotokolls).

Der Untergrund entsteht durch ein zweites, unabhaengig waehrend der Suchzeit
T_S einfallendes Myon, das einen falschen Stoppimpuls ausloest. Da ein solches
Myon zu jedem Zeitpunkt gleich wahrscheinlich eintrifft, ist der Untergrund
FLACH (konstant ueber alle Kanaele).

Er wird auf zwei Wegen bestimmt:
  Methode 1: y-Achsenabschnitt U aus dem Exponential-Fit (siehe lifetime.py).
  Methode 2: aus der Poisson-Wahrscheinlichkeit fuer ein zweites Myon.
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
# EINSTELLUNGEN  -- die mit (!) MUSST du an eure Messung anpassen
# ----------------------------------------------------------------------
SPE_PFAD   = "raw/2026_0611.Spe"

CAL_A = 0.021685       # us / Kanal   (aus channel.py)
CAL_B = 0.155200       # us
SKIP_FIRST = 0         # erste Kanaele nicht in den Fit einbeziehen

N_START = 986743.0     # (!) Startimpulse vom Zaehler (nicht in der .Spe gespeichert!)
T_S_US  = 20         # (!) Suchzeit / TAC-Messbereich in us.
                       #     None -> automatisch = MCA-Zeitbereich aus der Kalibrierung.
                       #     Altprotokoll nutzte hier 20 us. Setze euren echten Wert!


# ----------------------------------------------------------------------
# .Spe EINLESEN
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

    live, real = (int(x) for x in sektion("$MEAS_TIM:")[0].split())
    ch0, ch1   = (int(x) for x in sektion("$DATA:")[0].split())
    n_ch       = ch1 - ch0 + 1
    di         = lines.index("$DATA:")
    counts     = np.array([int(x) for x in lines[di + 2: di + 2 + n_ch]], dtype=float)
    return live, real, counts


live, real, counts = lies_spe(SPE_PFAD)
t_mess = real                                    # Messdauer in s (Real Time)

kanal = np.arange(len(counts))
t_us  = CAL_A * kanal + CAL_B
letzter_kanal = int(np.max(np.where(counts > 0)))
mca_bereich   = CAL_A * (len(counts) - 1) + CAL_B   # us

if T_S_US is None:
    T_S_US = mca_bereich
T_S = T_S_US * 1e-6                               # in s


# ----------------------------------------------------------------------
# METHODE 1: U aus dem Exponential-Fit
# ----------------------------------------------------------------------
einbezogen = (kanal >= SKIP_FIRST) & (kanal <= letzter_kanal)

def zerfall(t, N0, lam, U):
    return N0 * np.exp(-lam * t) + U

popt, pcov = curve_fit(zerfall, t_us[einbezogen], counts[einbezogen],
                       p0=[counts.max(), 0.46, 0.3], maxfev=10000)
perr = np.sqrt(np.diag(pcov))
U_fit = ufloat(popt[2], perr[2])


# ----------------------------------------------------------------------
# METHODE 2: Poisson-Statistik
# ----------------------------------------------------------------------
N_start = ufloat(N_START, np.sqrt(N_START))      # Poisson-Fehler auf die Zaehlung
N_myon  = N_start / t_mess                        # mittlere Einfallsrate [1/s]
mu      = T_S * N_myon                             # erwartete Myonen in T_S
P1      = mu * unp.exp(-mu)                         # P(genau ein weiteres Myon)
N_fehl  = N_start * P1                              # Zahl falscher Stopps
U_pois  = N_fehl / letzter_kanal                   # auf die befuellten Kanaele verteilt


# ----------------------------------------------------------------------
# AUSGABE
# ----------------------------------------------------------------------
abweichung = abs(U_pois.n - U_fit.n) / U_pois.n * 100

print("----- Eingangsgroessen -----")
print(f"t_Mess         = {t_mess} s")
print(f"N_start        = {N_start:.0f}   (Zaehleranzeige, NICHT aus .Spe!)")
print(f"T_S            = {T_S_US:.2f} us  (MCA-Bereich = {mca_bereich:.2f} us)")
print(f"befuellte Kanaele bis Kanal {letzter_kanal}\n")

print("----- Methode 1: Fit -----")
print(f"U_Fit          = {U_fit} Counts/Kanal\n")

print("----- Methode 2: Poisson -----")
print(f"N_Myon         = {N_myon} 1/s")
print(f"mu = T_S*N_Myon= {mu}")
print(f"P(1)           = {P1}")
print(f"N_fehl         = {N_fehl}")
print(f"U_Poisson      = {U_pois} Counts/Kanal\n")

print(f"Abweichung der beiden Methoden: {abweichung:.1f} %")