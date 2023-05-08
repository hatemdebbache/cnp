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
import MatVec

# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo1()
# ==============================================================
def exo1():
    U = np.zeros((4,4))
    a = 9000 ; b = 3000 ; c = 2000; v = 4000
    MatVec.remplir_matrice(0, 0, 1, 3, U, a)
    MatVec.remplir_matrice(1, 2, 0, 0, U, b)
    MatVec.remplir_matrice(3, 3, 1, 3, U, c)
    MatVec.remplir_matrice(1, 2, 1, 3, U, v)
    print(U)
    
    R = 1.
    kmax = 1
    e = 0.
    R2 = R*R
    c = 0.5/(R2 + 1)
    k = 0
    while (True):
        k = k + 1
        if (k > kmax):
            k = -1
            break
        
        propre = True
        for i in range(1,2+1):
            for j in range(1,2+1):        
                anc = U[i,j]
                U[i,j] = c*(U[i,j+1] + U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
                if (abs(anc - U[i,j]) > e):
                    propre = False
            
            j = 3
            anc = U[i,j]
            U[i,j] = c*(2*U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
            if (abs(anc - U[i,j]) > e):
                    propre = False

        if (propre == True):
            break
    
    print(k)
    print(U)
# ==============================================================
def exo2():
    pass
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
def f(x):
    pass
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