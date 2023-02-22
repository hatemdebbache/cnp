# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:14:25 2020

@author: user
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
plt.style.use('dark_background')
# plt.style.use('classic')
# ==============================================================
def FixeEchelle(xmin, xmax, ymin, ymax):
    plt.clf()
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

# ==============================================================
def TraceAxes(xo = 0.0, yo = 0.0, coulaxes='black', coulgr ='grey'):
    plt.grid(color = coulgr, linewidth = 0.5)
    plt.axvline(xo, color = coulaxes, linewidth = 1.0)
    plt.axhline(yo, color = coulaxes, linewidth = 1.0)

# ==============================================================
def TracePoints(x, y, relie=False, couleur = 'red', epaisseur = 1.0):
    if (relie):
        plt.plot(x, y, color = couleur, linewidth = epaisseur)
    else:
        plt.scatter(x, y, color = couleur, linewidths = epaisseur)

# ==============================================================
def TraceFonc(fonc, xmin, xmax, npts=101, relie=True, couleur='blue', epaisseur=1.0):
    x = [xmin + (xmax - xmin) * (k /(npts-1)) for k in range(npts)]
    y = [fonc(a) for a in x]
    TracePoints(x, y, relie=relie, couleur=couleur, epaisseur=epaisseur)

# ==============================================================
