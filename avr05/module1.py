# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020
Updated on Sat Feb 02 11:37:45 2023

@author: user
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi, cos, sin, log
# ----- User-Defined-Modules Importing
sys.path.insert(0, "./../BIBLIO/")
import cbase
import graphe
import edo
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
gALPHA = 1.
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo2()
# ==============================================================
def exo1():
    # Réglage du graphe
    ymin, ymax = -2., 4.
    zmin, zmax = -1., 4.
    graphe.FixeEchelle(ymin, ymax, zmin, zmax)
    graphe.TraceAxes()
    
    # Paramètres pour RK
    Nsub = 100 ; m = 1
    a = 0. ; b = 3.
    
    # Cercle
    r = 0.1 ; ntraj = 40 ; q = 2.*pi/ntraj
    Xc = [r*cos(k*q) for k in range(ntraj)]
    Yc = [r*sin(k*q) for k in range(ntraj)]
    
    # Les droite yp = 0 et zp = 0
    graphe.TraceFonc(f1, ymin, ymax, couleur='white', epaisseur=2.)
    plt.axvline(0., lw=2., c='white') # Droite vertical z = 0
    graphe.TraceFonc(f3, ymin, ymax, couleur='white', epaisseur=2.)
    graphe.TraceFonc(f4, ymin, ymax, couleur='white', epaisseur=2.)
    
    # PC1 : [0.,0.]   Source
    for xc, yc in zip(Xc, Yc):
        ya = xc ; za = yc
        c, t, y, z = edo.RungeKutta2(eq_diff, a, b, ya, za, Nsub, m)
        graphe.TracePoints(y, z, couleur='red', relie=True)

    # PC2 : [1., 2.]  Puit
    for xc, yc in zip(Xc, Yc):
        ya = 1. + xc ; za = 2. + yc
        c, t, y, z = edo.RungeKutta2(eq_diff, b, a, ya, za, Nsub, m)
        graphe.TracePoints(y, z, couleur='green', relie=True)

    # PC3 : [3., 0.]  Selle
    for xc, yc in zip(Xc, Yc):
        ya = 3. + xc ; za = yc
        c, t, y, z = edo.RungeKutta2(eq_diff, a, b, ya, za, Nsub, m)
        graphe.TracePoints(y, z, couleur='blue', relie=True)
        c, t, y, z = edo.RungeKutta2(eq_diff, b, a, ya, za, Nsub, m)
        graphe.TracePoints(y, z, couleur='blue', relie=True)

    # PC4 : [0., 1.]  Selle
    for xc, yc in zip(Xc, Yc):
        ya = xc ; za = 1. + yc
        c, t, y, z = edo.RungeKutta2(eq_diff, a, b, ya, za, Nsub, m)
        graphe.TracePoints(y, z, couleur='yellow', relie=True)
        c, t, y, z = edo.RungeKutta2(eq_diff, b, a, ya, za, Nsub, m)
        graphe.TracePoints(y, z, couleur='yellow', relie=True)
        
    plt.show()
# ==============================================================
def exo2():
    graphe.FixeEchelle(-2., 2., -2., 2.)
    graphe.TraceAxes()
    r = 0.1
    ntraj = 20
    q = 2*pi/ntraj
    for k in range(ntraj):
        ya = r*cos(k*q) ; za = 1. + r*sin(k*q)
        c, t, y, z = edo.RungeKutta2(samia, 0., 3., ya, za, 100, 100)
        graphe.TracePoints(y, z)
    for k in range(ntraj):
        ya = r*cos(k*q) ; za = -1. + r*sin(k*q)
        c, t, y, z = edo.RungeKutta2(samia, 0., 3., ya, za, 100, 100)
        graphe.TracePoints(y, z)
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 5., 1., 1., 100, 100)
    graphe.TracePoints(y, z, couleur='blue')
    
    ntraj = 10
    for k in range(2*ntraj):
        ya = 0. ; za = 1. - k/ntraj
        c, t, y, z = edo.RungeKutta2(samia, 0., 5., ya, za, 100, 100)
        graphe.TracePoints(y, z, couleur='cyan', relie=True)
    
    plt.show()
# ==============================================================
def exo3():
    pass
# ==============================================================
def exo4():
    pass
# ==============================================================
def exo5():
    pass
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def samia(t, y, z):
    return gALPHA - z*z , y
# ==============================================================
def eq_diff(t, y, z):
    return y*(3. - y - z), z*(1 + y - z) 
# ==============================================================
def f1(y):
    return 3. - y
# ==============================================================
def f2(z):
    return 0.
# ==============================================================
def f3(y):
    return 1. + y
# ==============================================================
def f4(y):
    return 0.
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================