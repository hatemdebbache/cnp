# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 13:15:34 2020
Updated on Sat Feb 02 11:37:45 2023

@author: user
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
#       selle       puit        source      spirale attractive
#                       EXERCICES
#
# ==============================================================
def main():
    return exo2()
# ==============================================================
def exo():
    # DONT EXECUTE
    n = 500 ; f = 0.01 
    # win = [-0.5*f*n, 0.5*f*n]
    # plt.xlim(win) ; plt.ylim(win)
    X = [] ; Y = []
    for i in range(n):
        for j in range(n):
            X.append(f*(i-n//2))
            Y.append(f*(j-n//2))
    a = 0. ; b = 5.
    Nsub = 1000 ; m = 10
    update = True
    if update:    
        with open('data.npy', 'wb') as file:
            for ya, za in zip(X, Y):
                c, t, y, z = edo.RungeKutta2(samia, a, b, ya, za, Nsub, m)
                # plt.plot(y, z)
                np.save(file, [y, z])
    else:
        from IPython import get_ipython # Used to access cache to reduce runtime
        ip = get_ipython()
        res = ip.user_ns['res']
        for instance in res:
            y, z = instance
            plt.plot(y, z)
    # plt.grid()
    # plt.show()    
# ==============================================================
def exo2():
    fig, ax = plt.subplots()
    ax.set_xlim((-2.,4.))
    ax.set_ylim((-1.,4.))
    Nsub = 100 ; m = 1
    a = 0. ; b = 3.
    # Cercle
    r = 0.1 ; npts = 40 ; theta = 2.*pi/npts
    Xc = [r*cos(k*theta) for k in range(npts)]
    Yc = [r*sin(k*theta) for k in range(npts)]
    y, z = [], []
    
    # PC1 : [0.,0.]   Source
    for xc, yc in zip(Xc, Yc):
        ya = xc ; za = yc
        c, t, y1, z1 = edo.RungeKutta2(samia, a, b, ya, za, Nsub, m)
        y = np.concatenate((y, y1)) ; z = np.concatenate((z, z1))
    # PC2 : [1., 2.]  Puit
    for xc, yc in zip(Xc, Yc):
        ya = 1. + xc ; za = 2. + yc
        c, t, y2, z2 = edo.RungeKutta2(samia, b, a, ya, za, Nsub, m)
        y = np.concatenate((y, y2)) ; z = np.concatenate((z, z2))
    # PC3 : [3., 0.]  Selle
    for xc, yc in zip(Xc, Yc):
        ya = 3. + xc ; za = yc
        c, t, y3, z3 = edo.RungeKutta2(samia, a, b, ya, za, Nsub, m)
        y = np.concatenate((y, y3)) ; z = np.concatenate((z, z3))
        c, t, y3, z3 = edo.RungeKutta2(samia, b, a, ya, za, Nsub, m)
        y = np.concatenate((y, y3)) ; z = np.concatenate((z, z3))
    # PC4 : [0., 1.]  Selle
    for xc, yc in zip(Xc, Yc):
        ya = xc ; za = 1. + yc
        c, t, y4, z4 = edo.RungeKutta2(samia, a, b, ya, za, Nsub, m)
        y = np.concatenate((y, y4)) ; z = np.concatenate((z, z4))
        c, t, y4, z4 = edo.RungeKutta2(samia, b, a, ya, za, Nsub, m)
        y = np.concatenate((y, y4)) ; z = np.concatenate((z, z4))
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 3., 2., 0., Nsub, m)
    plt.plot(y, z, color='cyan')
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 3., 2., -0.01, Nsub, m)
    plt.plot(y, z, color='cyan')
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 3., 4., 0., Nsub, m)
    plt.plot(y, z, color='cyan')
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 3., 4., 0.01, Nsub, m)
    plt.plot(y, z, color='cyan')
    
    c, t, y, z = edo.RungeKutta2(samia, 0., 3., -0.5, 2., Nsub, m)
    plt.plot(y, z, color='yellow')
    
    return y, z
# ==============================================================
#
#                       ADDITIONAL FUNCTIONS
#
# ==============================================================
def samia(t, y, z):
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
    
    y, z = main()
    
# ==============================================================