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
    exo4()
# ==============================================================
def exo1():
    U = np.zeros((11,11))
    a = 9000 ; b = 3000 ; c = 2000
    MatVec.remplir_matrice(0, 0, 4, 10, U, a)
    MatVec.remplir_matrice(10, 10, 1, 10, U, c)
    MatVec.remplir_matrice(1, 5, 3, 3, U, b)
    MatVec.remplir_matrice(5, 5, 0, 2, U, b)
    MatVec.remplir_matrice(6, 9, 0, 0, U, b)
    MatVec.remplir_matrice(1, 5, 4, 10, U, 5500)
    MatVec.remplir_matrice(6, 9, 1, 10, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    R = 1.2
    kmax = 250
    e = 0.01
    R2 = R*R
    c = 0.5/(R2 + 1)
    k = 0
    while (True):
        k = k + 1
        if (k > kmax):
            k = -1
            break
        
        propre = True
        for i in range(1,5+1):
            for j in range(4,9+1):        
                anc = U[i,j]
                U[i,j] = c*(U[i,j+1] + U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
                if (abs(anc - U[i,j]) > e):
                    propre = False
            j = 10
            anc = U[i,j]
            U[i,j] = c*(2*U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
            if (abs(anc - U[i,j]) > e):
                    propre = False

        for i in range(6,9+1):
            for j in range(1,9+1):
                anc = U[i,j]
                U[i,j] = c*(U[i,j+1] + U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
                if (abs(anc - U[i,j]) > e):
                    propre = False
            j = 10
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
    U = np.zeros((11,11))
    a = 9000 ; b = 3000 ; c = 2000 ; d = 1000
    MatVec.remplir_matrice(0, 0, 1, 9, U, a)
    MatVec.remplir_matrice(10, 10, 1, 9, U, c)
    MatVec.remplir_matrice(1, 9, 0, 0, U, b)
    MatVec.remplir_matrice(1, 9, 10, 10, U, d)
    MatVec.remplir_matrice(1, 9, 1, 9, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    
    R = 1.2
    kmax = 250
    e = 0.001
    R2 = R*R
    c = 0.5/(R2 + 1)
    k = 0
    while (True):
        k = k + 1
        if (k > kmax):
            k = -1
            break
        
        propre = True
        for i in range(1,9+1):
            for j in range(1,9+1):        
                anc = U[i,j]
                U[i,j] = c*(U[i,j+1] + U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
                if (abs(anc - U[i,j]) > e):
                    propre = False

        if (propre == True):
            break
    
    print(k)
    print(U)
# ==============================================================
def exo3():
    U = np.zeros((6,6))
    a = 9000 ; b = 3000 ; c = 2000 ; d = 1000
    MatVec.remplir_matrice(0, 0, 1, 4, U, a)
    MatVec.remplir_matrice(5, 5, 1, 4, U, c)
    MatVec.remplir_matrice(1, 4, 0, 0, U, b)
    MatVec.remplir_matrice(1, 4, 5, 5, U, d)
    MatVec.remplir_matrice(1, 4, 1, 4, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    
    R = 1.2
    kmax = 250
    e = 0.01
    R2 = R*R
    c = 0.5/(R2 + 1)
    k = 0
    while (True):
        k = k + 1
        if (k > kmax):
            k = -1
            break
        
        propre = True
        for i in range(1,4+1):
            for j in range(1,4+1):        
                anc = U[i,j]
                U[i,j] = c*(U[i,j+1] + U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
                if (abs(anc - U[i,j]) > e):
                    propre = False

        if (propre == True):
            break
    
    print(k)
    print(U)
# ==============================================================
def exo4():
    U = np.zeros((11,11))
    a = 9000 ; b = 3000 ; c = 2000 ; d = 1000
    MatVec.remplir_matrice(0, 0, 1, 9, U, a)
    MatVec.remplir_matrice(10, 10, 1, 9, U, c)
    MatVec.remplir_matrice(1, 9, 0, 0, U, b)
    MatVec.remplir_matrice(1, 9, 10, 10, U, d)
    MatVec.remplir_matrice(1, 9, 1, 9, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    
    R = 1.2
    kmax = 250
    e = 0.001
    R2 = R*R
    c = 0.5/(R2 + 1)
    w = 1.55
    k = 0
    while (True):
        k = k + 1
        if (k > kmax):
            k = -1
            break
        
        propre = True
        for i in range(1,9+1):
            for j in range(1,9+1):        
                anc = U[i,j]
                nouv = c*(U[i,j+1] + U[i,j-1] + R2*(U[i-1,j] + U[i+1,j]))
                U[i,j] = anc + w*(nouv - anc)
                if (abs(anc - U[i,j]) > e):
                    propre = False

        if (propre == True):
            break
    
    print(k)
    print(U)
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