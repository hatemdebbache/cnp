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
import approx
import edo
import MatVec as mat
import cbase

# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
cs = 6         # Numbre of most significant digits to be printed
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo1()
    # exo2()
    # exo3()
    # exo4()
    # exo5()
# ==============================================================
def exo1():
    print("==================== EXO 1 =====================")
    # Données de l'exercice
    a = 8e3*(1 + gBETA) ; b = 1e3*(1 + gBETA) ; c = 5e3*(1 + gBETA)
    P = 9e3 ; v = 5e3
    R = 1.2
    # Creation de U
    u = np.zeros((11,11))
    # Valeurs initiales
    mat.remplir_matrice(1, 6, 0, 9, u, v)
    mat.remplir_matrice(7, 10, 5, 9, u, v)
    # Conditions de Dirichlet
    mat.remplir_matrice(0,0,0,9,u,a)
    mat.remplir_matrice(1,10,10,10,u,c)
    mat.remplir_matrice(7,7,0,4,u,b)
    mat.remplir_matrice(7,10,4,4,u,b)
    u[5,7] = P
    # Resolution
    R2 = R*R
    c = 0.5/(R2 + 2)
    for w in [1., 1.2, 1.5]:
        kmax = 3 ; e = 0.
        k = 0
        while (True):
            k += 1
            if (k>kmax):
                k = -1 ; break
            
            done = True
            # Premier block
            for i in range(1,6+1):
                # Premiere colonne
                anc = u[i,0]
                nouv = c*(4*u[i,1] + R2*(1-1.5/i)*u[i-1,0] + R2*(1+1.5/i)*u[i+1,0])
                u[i,0] = anc + w*(nouv-anc)
                if (abs(anc - u[i,0]) > e):
                    done = False
                # Le rest
                for j in range(1,9+1):
                    if (i==5 and j==7): 
                        continue  # Ne modifier pas P
                    anc = u[i,j]
                    nouv = c*(2*(u[i,j-1] + u[i,j+1]) + R2*(1-1.5/i)*u[i-1,j] + R2*(1+1.5/i)*u[i+1,j])
                    u[i,j] = anc + w*(nouv-anc)
                    if (abs(anc - u[i,j]) > e):
                        done = False
            # Deuxieme block
            for i in range(7,9+1):
                for j in range(5,9+1):
                    anc = u[i,j]
                    nouv = c*(2*(u[i,j-1] + u[i,j+1]) + R2*(1-1.5/i)*u[i-1,j] + R2*(1+1.5/i)*u[i+1,j])
                    u[i,j] = anc + w*(nouv-anc)
                    if (abs(anc - u[i,j]) > e):
                        done = False
            # Dernier ligne
            for j in range(5,9+1):
                anc = u[10,j]
                nouv = c*(2*(u[10,j-1] + u[10,j+1]) + 2*R2*u[9,j])
                u[10,j] = anc + w*(nouv-anc)
                if (abs(anc - u[10,j]) > e):
                    done = False
            
            if (done):
                break
    
        V1 = u[5,0] ; V2 = u[10,7]
        print(f'w = {w} \tcout = {k}\nV1 = {V1:.0f}\tV2 = {V2:.0f}\n')
        
    # edp.heatmap(u)
    print("================================================\n")
# ==============================================================
def exo2():
    print("==================== EXO 2 =====================")
    xmin = -1 ; xmax = 0
    ymin = -1.5 ; ymax = 1.5
    N = 51 ; cs = 4
    x = np.linspace(xmin,xmax,N)
    y = np.zeros(N)
    for i in range(N):
        y[i] = f(x[i])
    
    a, b = approx.Droite_MC(x, y)
    y_approch = a*x + b
    print(f'a = {a:.{cs}g}\tb = {b:.{cs}g}\n')
    S = approx.EQM(x, y, Y_approch=y_approch)
    print(f'S = {S:.{cs}g}\n')
    
    # graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    # graphe.TraceAxes()
    # graphe.TraceFonc(f, xmin, xmax)
    # graphe.TracePoints(x, y)
    # plt.show()
    print("================================================\n")
# ==============================================================
def exo3():
    print("==================== EXO 3 =====================")
    a = 1.   ; b = 2.
    ya = 1.5 ; yb = 2.
    
    Nsub = 100
    x, y = edo.DFC2_GT(PQR, a, b, ya, yb, Nsub, True, True)
    imin = np.argmin(y)
    print(f'xmin = {x[imin]:.2g}\tymin= {y[imin]:.2g}\n')
    
    x1 = 1.3 ; x2 = 1.295
    i1 = 0 ; i2 = 0
    while (x[i1] < x1):
        i1 += 1
    while(x[i2] < x2):
        i2 += 1
    
    w = y[i1]
    print(f'w = {w:.{cs}g}\n')
    
    A1 = cbase.Poly_Lagrange(x, y, x2)
    print(f'A1 = {A1:.{cs}g}\n')
    
    
    Nsub = 400
    x, y = edo.DFC2_GT(PQR, a, b, ya, yb, Nsub, True, True)
    
    x2 = 1.295 ; i2 = 0
    while(x[i2] < x2):
        i2 += 1
    
    A2 = cbase.Poly_Lagrange(x, y, x2)
    print(f'A2 = {A2:.{cs}g}\n')
    
    print(f'E = {abs(A2-A1):.2g}\n')
    # xmin = 1 ; xmax = 2
    # ymin = 1 ; ymax = 2
    # graphe.FixeEchelle(xmin, xmax, ymin, ymax)
    # graphe.TraceAxes()
    # graphe.TracePoints(x,y)
    # plt.show()
    print("================================================\n")
# ==============================================================
def exo4():
    print("==================== EXO 4 =====================")
    
    
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
def f(x):
    x = 0.5*pi*x + gBETA
    return cos(x) + sin(x)
# ==============================================================
def PQR(x):
    c = 1/(1 + gBETA*x*x)
    p = c
    q = c*x
    r = c*exp(-x)
    return p, q, r
# ==============================================================
#
#                       EXECUTION BLOCK
#
# ==============================================================
if (__name__ == "__main__"):
    
    main()
    
# ==============================================================