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
import syslin
import MatVec
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
# ==============================================================
#
#                       EXERCICES
#
# ==============================================================
def main():
    exo3()
# ==============================================================
def exo1():
    A = np.array([[4., -1., -1., 0.],
                  [-1., 4., 0., -1.],
                  [-1., 0., 4., -1.],
                  [0., -1., -1., 4.]])
    b = np.array([12000., 14000., 5000., 7000.])
    x = syslin.Elim_Gauss0(A, b)
    print(x)
# ==============================================================
def exo2():
    U = np.zeros((4,4))
    MatVec.remplir_matrice(0, 0, 1, 2, U, 9000)
    MatVec.remplir_matrice(3, 3, 1, 2, U, 2000)
    MatVec.remplir_matrice(1, 2, 0, 0, U, 3000)
    MatVec.remplir_matrice(1, 2, 3, 3, U, 5000)
    MatVec.remplir_matrice(1, 2, 1, 2, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    for k in range(5):
        for i in range(1,2+1):
            for j in range(1,2+1):
                U[i,j] = 0.25*(U[i,j+1] + U[i,j-1] + U[i-1,j] + U[i+1,j])
    print(U)
# ==============================================================
def exo3():
    U = np.zeros((11,11))
    MatVec.remplir_matrice(0, 0, 1, 9, U, 9000)
    MatVec.remplir_matrice(10, 10, 1, 9, U, 2000)
    MatVec.remplir_matrice(1, 9, 0, 0, U, 3000)
    MatVec.remplir_matrice(1, 9, 10, 10, U, 5000)
    MatVec.remplir_matrice(1, 9, 1, 9, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    for k in range(100):
        for i in range(1,9+1):
            for j in range(1,9+1):
                U[i,j] = 0.25*(U[i,j+1] + U[i,j-1] + U[i-1,j] + U[i+1,j])
    print(U)
# ==============================================================

# 4 & 5 added on apr26th

def exo4():
    U = np.zeros((11,11))
    MatVec.remplir_matrice(0, 0, 5, 9, U, 9000)
    MatVec.remplir_matrice(10, 10, 1, 9, U, 2000)
    MatVec.remplir_matrice(4, 9, 0, 0, U, 3000)
    MatVec.remplir_matrice(1, 4, 4, 4, U, 3000)
    MatVec.remplir_matrice(4, 4, 1, 4, U, 3000)
    MatVec.remplir_matrice(1, 9, 10, 10, U, 5000)
    MatVec.remplir_matrice(1, 9, 5, 9, U, 5500)
    MatVec.remplir_matrice(5, 9, 1, 4, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    for k in range(100):
        for i in range(1,4+1):
            for j in range(5,9+1):
                U[i,j] = 0.25*(U[i,j+1] + U[i,j-1] + U[i-1,j] + U[i+1,j])
        
        for i in range(5,9+1):
            for j in range(1,9+1):
                U[i,j] = 0.25*(U[i,j+1] + U[i,j-1] + U[i-1,j] + U[i+1,j])
    print(U)
# ==============================================================
def exo5():
    U = np.zeros((11,11))
    MatVec.remplir_matrice(0, 0, 5, 10, U, 9000)
    MatVec.remplir_matrice(10, 10, 1, 10, U, 2000)
    MatVec.remplir_matrice(4, 9, 0, 0, U, 3000)
    MatVec.remplir_matrice(1, 4, 4, 4, U, 3000)
    MatVec.remplir_matrice(4, 4, 1, 4, U, 3000)
    MatVec.remplir_matrice(1, 9, 10, 10, U, 5000)
    MatVec.remplir_matrice(1, 9, 5, 10, U, 5500)
    MatVec.remplir_matrice(5, 9, 1, 4, U, 5500)
    np.set_printoptions(precision=0)
    # print(U)
    for k in range(100):
        for i in range(1,4+1):
            for j in range(5,9+1):
                U[i,j] = 0.25*(U[i,j+1] + U[i,j-1] + U[i-1,j] + U[i+1,j])
            j = 10
            U[i,j] = 0.25*(2*U[i,j-1] + U[i-1,j] + U[i+1,j])
        
        for i in range(5,9+1):
            for j in range(1,9+1):
                U[i,j] = 0.25*(U[i,j+1] + U[i,j-1] + U[i-1,j] + U[i+1,j])
            j = 10
            U[i,j] = 0.25*(2*U[i,j-1] + U[i-1,j] + U[i+1,j])
            
    print(U)
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