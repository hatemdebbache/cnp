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
import graphe
import edo
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo2()
# ==============================================================
def exo1():
    to = 0. ; tf = 10.
    
    Nsub = 1000 ; m = 10
    ntraj = 10  ; r = 0.1
    ymin, ymax = -5., 5.
    zmin, zmax = -5., 5.
    graphe.FixeEchelle(ymin, ymax, zmin, zmax)
    graphe.TraceAxes()
    PTC = [[0.,0.], [0.,1.], [-1.,-2.]]
    for PC in PTC:
        draw_trajectoire_pc(PC, yasmine, to, tf, Nsub, m, r, ntraj)
    
    to, tf = tf, to
    for PC in PTC:
        draw_trajectoire_pc(PC, yasmine, to, tf, Nsub, m, r, ntraj)
    
    plt.show()
# ==============================================================
def exo2():
    to = 0. ; tf = 5.
    
    Nsub = 1000 ; m = 10
    ntraj = 20  ; r = 0.1
    ymin, ymax = -5., 5.
    zmin, zmax = -5., 5.
    graphe.FixeEchelle(ymin, ymax, zmin, zmax)
    graphe.TraceAxes()
    PTC = [[0.,0.], [0.,1.], [0.5,0.5]]
    for PC in PTC:
        draw_trajectoire_pc(PC, eq_diff, to, tf, Nsub, m, r, ntraj)
    
    to, tf = tf, to
    for PC in PTC:
        draw_trajectoire_pc(PC, eq_diff, to, tf, Nsub, m, r, ntraj)
    
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
def draw_trajectoire_pc(PC, eq_diff, to, tf, Nsub, m, r, ntraj):
    Y, Z = np.empty(0), np.empty(0)
    Yc, Zc = PC
    q = 2*pi/ntraj    
    for k in range(ntraj):
        ya = Yc + r*cos(k*q)
        za = Zc + r*sin(k*q)
        c, t, y, z = edo.RungeKutta2(eq_diff, to, tf, ya, za, Nsub, m)
        graphe.TracePoints(y, z, relie=True)
# ==============================================================
def yasmine(t, y, z):
    return z*(1. - y*y), y*(1. - y + z)
# ==============================================================
def eq_diff(t, y, z):
    return y*(1. - y - z), z*(1. - 2.*y)
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================