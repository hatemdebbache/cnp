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
import cbase
import graphe
import edo

# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
cs = 9         # Numbre of most significant digits to be printed
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    # exo1()
    # exo2()
    exo3()
    # exo4()
    # exo5()
# ==============================================================
def exo1():
    print("==================== EXO 1 =====================")
    # Données du problème
    yo = 10*gBETA
    a = 0. ; b = 2.
    # ------ Partie 1
    H = 0.1 ; h = 5e-2
    Nsub = int((b - a)/H) ; m = int(H/h)
    
    # Solution exacte
    V = exacte(1.3)
    # Resolution numérique
    c, t, VE = edo.Euler1(eq_diff_1, a, b, yo, Nsub, m)
    c, t, VP = edo.PointMilieu1(eq_diff_1, a, b, yo, Nsub, m)
    c, t, VR = edo.RungeKutta1(eq_diff_1, a, b, yo, Nsub, m)
    # indice de l'instant t = 1.3s
    i = int((1.3 - a)/H)
    EE = abs(VE - V) ; EP = abs(VP - V) ; ER = abs(VR - V)
    print(f'Pour t = {t[i]}s:\nV = {V:.{cs}g}')
    # Par la méthode d'Euler
    print(f'Ve = {VE[i]:.{cs}g} \tEe = {EE[i]:.{cs}g}')
    # Par la méthode du Point Milieu
    print(f'Vpm = {VP[i]:.{cs}g} \tEpm = {EP[i]:.{cs}g}')
    # Par la méthode de Runge-Kutta
    print(f'Vrk = {VR[i]:.{cs}g} \tEr = {ER[i]:.{cs}g}\n')
    
    # ------ Pour 0.5*hi
    h *= 0.5
    Nsub = int((b - a)/H) ; m = int(H/h)
    c, t, VE = edo.Euler1(eq_diff_1, a, b, yo, Nsub, m)
    c, t, VR = edo.RungeKutta1(eq_diff_1, a, b, yo, Nsub, m)
    i = int((1.3 - a)/H)
    EEn = abs(VE - V) ; ERn = abs(VR - V)
    # Par la méthode d'Euler
    print(f'Ve = {VE[i]:.{cs}g} \tEe = {EEn[i]:.{cs}g}')
    # Par la méthode de Runge-Kutta
    print(f'Vrk = {VR[i]:.{cs}g} \tEr = {ERn[i]:.{cs}g}\n')
    
    p1 = EEn/EE ; p2 = ERn/ER
    ne = abs(np.log2(p1)) ; nr = abs(np.log2(p2))
    print(f'ne = {ne[i]:.{cs}g} \tnr = {nr[i]:.{cs}g}\n')
    print("================================================\n")
# ==============================================================
def exo2():
    print("==================== EXO 2 =====================")
    # Données du problème
    L = gBETA ; C = 1e-2 ; R = 10. ; E = 5.
    a = 0. ; b = 0.7
    uo = 0. ; vo = 0.
    dt = 1e-4
    Nsub = 700 ; m = 10 # int((b-a)/(dt*Nsub))
    
    # Système d'équations différentielles
    def eq_diff_2(t, u, v):
        up = v
        vp = -v*R/L -u/(L*C) + E/(L*C)
        return up, vp
    c, t, u, v = edo.RungeKutta2(eq_diff_2, a, b, uo, vo, Nsub, m)
    
    t6 = 0.6 ; i6 = int((t6 - a)/(dt*m))
    usat = u[i6]
    print(f'tsat = {t[i6]:.{cs}g}s \tUsat = {usat:.{cs}g}V\n')
    
    i95 = 0
    while(u[i95] < 0.95*usat):
        i95 += 1
    csn = 4
    print(f't95 = {t[i95]:.{csn}g}s \tu95 = {u[i95]:.{csn}g}V\n')
    
    t95 = 3.*R*C ; i95 = int((t95 - a)/(dt*m))
    print(f't95 = {t[i95]:.{csn}g}s \t\tu95 = {u[i95]:.{csn}g}V\n')
    
    def eq_diff_3(t, q, i):
        qp = i
        ip = -(q/C + R*i)/L
        return qp, ip
    a = 0. ; b = 1.
    qo = 0. ; io = 30.
    m = 1 ; Nsub = int((b - a)/dt)
    
    c, t, q, i = edo.RungeKutta2(eq_diff_3, a, b, qo, io, Nsub, m)
    i4 = int((0.4 - a)/(dt*m))
    print(f't0.4 = {t[i4]:.{cs}g}s \tI0.4 = {i[i4]:.{cs}g}A\n')
    
    # graphe.FixeEchelle(a, b, 0., 6.)
    # graphe.TraceAxes()
    # graphe.TracePoints(t, u)
    # plt.show()
    print("================================================\n")
# ==============================================================
def exo3():
    print("==================== EXO 3 =====================")
    csn = 4
    # Données
    a = 1. ; b = 1. ; c = gBETA ; d = 1.
    dt = 0.1
    yo = 5. ; zo = 1.
    a = 0. ; b = 100.
    m = 1 ; Nsub = int((b - a)/dt)
    
    # Eq Diff
    # def Lotka(t, y, z):
    #     yp = a*y*(1 - z/b)
    #     zp = -c*z*(1 - y/d)
    #     return yp, zp
    def Lotka(t, y, z):
        return y*(1 - z) , -gBETA*z*(1 - y)
    c, t, y, z = edo.PointMilieu2(Lotka, a, b, yo, zo, Nsub, m)
    
    a3, b3 = 400, 500
    a4, b4 = 600, 700
    y13 = max(y[a3:b3]) ; y14 = max(y[a4:b4])
    z13 = max(z[a3:b3]) ; z14 = max(z[a4:b4])
    print(f'3ème pique:\ny1 = {y13:.{csn}g} \ty2 = {z13:.{csn}g}\n')
    print(f'4ème pique:\ny1 = {y14:.{csn}g} \ty2 = {z14:.{csn}g}\n')
    
    hmin, hmax = np.argmin(y), np.argmax(y)
    vmin, vmax = np.argmin(z), np.argmax(z)
    
    # print(f'y1 = {y[hmin]:.{csn}g} \t{y[hmax]:.{csn}g}')
    # print(f'y2 = {z[vmin]:.{csn}g} \t{z[vmax]:.{csn}g}\n')
    
    print(f'y1 = {y[vmin]:.{csn}g} \t{y[vmax]:.{csn}g}')
    print(f'y2 = {z[hmin]:.{csn}g} \t{z[hmax]:.{csn}g}\n')
    
    # graphe.FixeEchelle(a, b, 0., 6.)
    # graphe.TracePoints(t, y, couleur='blue')
    # graphe.TracePoints(t, z, couleur='red')
    # graphe.FixeEchelle(0., 6., 0., 3.)
    # graphe.TraceAxes()
    # graphe.TracePoints(y, z)
    # plt.show()
    
    print("================================================\n")
# ==============================================================
def exo4():
    print("==================== EXO 4 =====================")
    # Données 
    L = 1 + gBETA ; g = 9.81
    theto = 2*gBETA ; omo = 0.
    a = 0. ; b = 10.
    dt = 0.05
    Nsub = 5 ; m = int((b-a)/(dt*Nsub))
    # m = 1 ; Nsub = int((b - a)/dt)
    print(m)
    def eq_diff_pondule(t, thet, om):
        thetp = om
        omp = -g*sin(thet)/L
        return thetp, omp
    c, t, thet, om = edo.RungeKutta2(eq_diff_pondule, a, b, theto, omo, Nsub, m)
    print("Tableau de theta:\n")
    cbase.affiche_tableau(np.array([t, thet]).T, cs)
    print()
    omega = sqrt(g/L)
    
    def sol_approch(t):
        return 2*gBETA*cos(omega*t)
    
    graphe.FixeEchelle(a, b, 0., 0.3)
    graphe.TraceAxes()
    graphe.TracePoints(t, thet)
    graphe.TraceFonc(sol_approch, a, b)
    plt.show()
    print("================================================\n")
# ==============================================================
def exo5():
    print("==================== EXO 5 =====================")
    
    
    print("================================================\n")
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def eq_diff_1(t, y):
    return -2.*y
# ==============================================================
def exacte(t):
    return 10*gBETA*exp(-2.*t)
# ==============================================================

# ==============================================================
def g(x):
    pass
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================