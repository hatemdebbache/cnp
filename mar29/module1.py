# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020
Updated on Sat Feb 02 11:37:45 2023

@author: user
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi, cos, sin, log, sqrt
# ----- User-Defined-Modules Importing
sys.path.insert(0, "./../BIBLIO/")
import graphe
import edo
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
gALPHA = 0.5
gOMEGA = 1.
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo4()
# =============================================================
def exo1():
    graphe.FixeEchelle(-5., 5., -5., 5.)
    graphe.TraceAxes()
    graphe.TraceFonc(samir, -5., -1.001, epaisseur=3., couleur='red', npts=1000)
    graphe.TraceFonc(samir, -0.999, 0.999, epaisseur=3., couleur='red', npts=1000)
    graphe.TraceFonc(samir, 1.001, 5., epaisseur=3., couleur='red', npts=1000)
    graphe.TraceFonc(hatem, -5., 5., couleur='blue', epaisseur=3.)
    
    c, t, y, z = edo.RungeKutta2(van_der_pol, 0., 20., -3., 3., 1000, 100)
    graphe.TracePoints(y, z, epaisseur=2., couleur='green')
    
    c, t, y, z = edo.RungeKutta2(van_der_pol, 0., 20., -4., 3., 1000, 100)
    graphe.TracePoints(y, z, epaisseur=2., couleur='green')
    
    c, t, y, z = edo.RungeKutta2(van_der_pol, 0., 20., 3., -3., 1000, 100)
    graphe.TracePoints(y, z, epaisseur=2., couleur='yellow')
    
    c, t, y, z = edo.RungeKutta2(van_der_pol, 0., 20., 1.8, -2., 1000, 100)
    graphe.TracePoints(y, z, epaisseur=2., couleur='yellow')
    
    c, t, y, z = edo.RungeKutta2(van_der_pol, 0., 20., 0., 4., 1000, 100)
    graphe.TracePoints(y, z, epaisseur=2., couleur='cyan')
    
    c, t, y, z = edo.RungeKutta2(van_der_pol, 0., 20., 0.9, 0.9, 1000, 100)
    graphe.TracePoints(y, z, epaisseur=2., couleur='magenta')
    
    plt.show()
# ==============================================================
def exo2():
    graphe.FixeEchelle(0., 20., -5., 5.)
    graphe.TraceAxes()
    
    c, t, y, z = edo.RungeKutta2(van_der_pol, 0., 20., 0.9, 0.9, 1000, 100)
    graphe.TracePoints(t, y, epaisseur=2., couleur='magenta')
    
    plt.show()
# ==============================================================
def exo3():
    graphe.FixeEchelle(-2., 2., -2., 2.)
    graphe.TraceAxes()
    
    c, t, r, q = edo.RungeKutta2(samia, 0., 20., 0.1, 0., 400, 100)
    y = np.zeros(401)
    z = np.zeros(401)
    for k in range(0, 401):
        y[k] = r[k]*cos(q[k])
        z[k] = r[k]*sin(q[k])
    graphe.TracePoints(y, z)
    
    c, t, r, q = edo.RungeKutta2(samia, 0., 20., 1.9, 0., 400, 100)
    y = np.zeros(401)
    z = np.zeros(401)
    for k in range(0, 401):
        y[k] = r[k]*cos(q[k])
        z[k] = r[k]*sin(q[k])
    graphe.TracePoints(y, z, couleur='blue')
    
    plt.show()
# ==============================================================
def exo4():
    graphe.FixeEchelle(-2., 2., -2., 2.)
    graphe.TraceAxes()
    
    for k in range(30):
        qi = 2*k*pi/30
        c, t, r, q = edo.RungeKutta2(samia2, 0., 6., 1.8, qi, 400, 100)
        y = np.zeros(401)
        z = np.zeros(401)
        for k in range(0, 401):
            y[k] = r[k]*cos(q[k])
            z[k] = r[k]*sin(q[k])
        graphe.TracePoints(y, z)
        
    plt.show()
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def van_der_pol(t, y, z):
    return z, -y - gALPHA*(y*y - 1.)*z 
# ==============================================================
def samir(y):
    q = 1. - y*y
    if not (q == 0):    
        return y/(gALPHA*q)
    else:
        return 0.
# ==============================================================
def hatem(x):
    return 0.
# ==============================================================
def samia(t, r, q):
    return r*(1. - r), gALPHA
# ==============================================================
def samia2(t, r, q):
    return r*(1. - r), gALPHA - sin(q)
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================