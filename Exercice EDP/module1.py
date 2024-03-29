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
import MatVec
import disque
import ctrl
import graphe

# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
cs = 4
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo2()
# ==============================================================
def exo1():
    # Parametres de l'exercice
    R = 1
    a = 1200*sqrt(gBETA) ; b = 200/gBETA ; p = 1500 ; c = 1000
    v_init = 1000
    kmax = 1000 ; e = 1e-4 ; w = 1.2
    
    # Remplissage du matrice initiale
    u = np.zeros((24,15))
    MatVec.remplir_matrice(0, 0, 0, 14, u, a)
    MatVec.remplir_matrice(1, 9, 0, 0, u, b)
    MatVec.remplir_matrice(1, 9, 14, 14, u, c)
    MatVec.remplir_matrice(1, 9, 1, 13, u, v_init)
    MatVec.remplir_matrice(10, 23, 0, 0, u, b)
    u[5,7] = p

    for i in range(10,24):
        u[i,-(i-10)+14] = b
        for j in range(1,24-i):
            u[i,j] = v_init
    
    # Resolution de l'EDP
    R2 = 2*R*R ; R14 = 0.25*R
    c = 0.5/(R2 + 1)
    k = 0
    while (True):
        # Verification de cout
        k = k + 1
        if (k > kmax):
            k = -1
            break
        done = True
        
        # 1er noyau
        for i in range(1,9+1):
            for j in range(1,13+1):
                # Le point (5,7) doit rester inchange
                if (i==5 and j ==7):
                    continue

                anc = u[i,j]
                nouv = c*(R14*(u[i-1,j-1] - u[i-1,j+1] - u[i+1,j-1] + u[i+1,j+1]) + R2*(u[i-1,j] + u[i+1,j]) + u[i,j-1] + u[i,j+1])
                u[i,j] = anc + w*(nouv - anc)
                if (abs(anc - u[i,j]) > e):
                    done = False 
        
        # 2eme noyau
        for i in range(10,23+1):
            for j in range(1,24-i):
                anc = u[i,j]
                nouv = c*(R2*(u[i-1,j] + u[i+1,j]) + u[i,j-1] + u[i,j+1])
                u[i,j] = anc + w*(nouv - anc)
                if (abs(anc - u[i,j]) > e):
                    done = False
        
        # Verification de la convergence
        if (done):
            break
    
    # Affichage
    print(f'U(5,6) = {u[5,6]:.{cs}g}\tU(9,14) = {u[9,14]:.{cs}g}\tU(18,4) = {u[18,4]:.{cs}g}')
    print(f'Cout = {k}')
    
    disque.EcrireMat(u, 'zaki.txt', 4)

# ==============================================================
def exo2():
    # La methode la plus simple de trouver le cout optimal est la methode graphique.
    # On trace k en fonction de w et trouver le minimum (NB: w appartient à [1,2[).
    # On peut retrouver la valeur du minimum en utilisant les modules 'derive' et 'racine'.
    # Le minimum correspond à la racine de d(resolution_edp)/dw = 0.
    # Il est preferable de localiser ce minimum graphiquemenet (le graphe avec faible resolution),
    #   puis utiliser la deuxieme dans le domain restreint pour obtenir plus de chiffres significatifs
    
    # graphe.FixeEchelle(0.9, 2.1, 0, 200)
    # graphe.TraceAxes()
    # graphe.TraceFonc(resolution_edp, 1.001, 1.999, npts=998)
    # plt.show()
    
    # D'apres le graphe
    print(f"Cout_min = 42 ; w_opt = 1.62 (3CS)")
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
def resolution_edp(w):
    # Parametres de l'exercice
    R = 1
    a = 1200*sqrt(gBETA) ; b = 200/gBETA ; p = 1500 ; c = 1000
    v_init = 1000
    kmax = 1000 ; e = 1e-4
    
    # Remplissage du matrice initiale
    u = np.zeros((24,15))
    MatVec.remplir_matrice(0, 0, 0, 14, u, a)
    MatVec.remplir_matrice(1, 9, 0, 0, u, b)
    MatVec.remplir_matrice(1, 9, 14, 14, u, c)
    MatVec.remplir_matrice(1, 9, 1, 13, u, v_init)
    MatVec.remplir_matrice(10, 23, 0, 0, u, b)
    u[5,7] = p

    for i in range(10,24):
        u[i,-(i-10)+14] = b
        for j in range(1,24-i):
            u[i,j] = v_init
    
    # Resolution de l'EDP
    R2 = 2*R*R ; R14 = 0.25*R
    c = 0.5/(R2 + 1)
    k = 0
    while (True):
        # Verification de cout
        k = k + 1
        if (k > kmax):
            k = -1
            break
        done = True
        
        # 1er noyau
        for i in range(1,9+1):
            for j in range(1,13+1):
                # Le point (5,7) doit rester inchange
                if (i==5 and j ==7):
                    continue

                anc = u[i,j]
                nouv = c*(R14*(u[i-1,j-1] - u[i-1,j+1] - u[i+1,j-1] + u[i+1,j+1]) + R2*(u[i-1,j] + u[i+1,j]) + u[i,j-1] + u[i,j+1])
                u[i,j] = anc + w*(nouv - anc)
                if (abs(anc - u[i,j]) > e):
                    done = False 
        
        # 2eme noyau
        for i in range(10,23+1):
            for j in range(1,24-i):
                anc = u[i,j]
                nouv = c*(R2*(u[i-1,j] + u[i+1,j]) + u[i,j-1] + u[i,j+1])
                u[i,j] = anc + w*(nouv - anc)
                if (abs(anc - u[i,j]) > e):
                    done = False
        
        # Verification de la convergence
        if (done):
            break
        
    return k
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