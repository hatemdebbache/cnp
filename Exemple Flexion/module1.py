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
import edo
# ----- Global Variables/Constants
gBETA = 0.1433 # PERSONAL CODE
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
    plt.clf()
    plt.grid()
    plt.xlabel('position (m)', fontsize=20.)
    plt.ylabel('Flèche (m)', fontsize=20.)
    
    # Input
    F = 2.        # N
    a = 1.35      # m
    L = 10.       # m
    EI =  2e8*(.2**4/6) # N.m² square cross-section of 20x20cm
    # Computing analytic sol
    x, y = simple_bend(F, a, L, EI)    
    # Plotting results
    plt.plot(x, y, linewidth=6., color='black')
    ls = '--'; clr='pink'; lw=4.
    plt.axvline(x=0, linestyle=ls, color=clr,linewidth=lw)
    plt.axvline(x=a, linestyle=ls, color=clr,linewidth=lw)
    plt.axvline(x=L, linestyle=ls, color=clr,linewidth=lw)

    # Computing solution with DFC2_GT
    def Moment(x):
        if (x < a):
            M = 1 - a/L
        else:
            M = a*(L-x)/L
        r = M*F/EI
        
        return 0., 0., r
    Nsub = 100
    x, y = edo.DFC2_GT(Moment, 0., L, 0., 0., Nsub, True, True)
    # Plotting as well
    plt.scatter(x, y, linewidth = 3., color='red')
    
    # Verif par flech max
    b = L - a
    if (a >= 0.5*L):
        ymax = -F*b*(3*L*L - 4*b*b)/(48*EI)
        plt.axhline(y=ymax, label='flèche-max',color='blue')
        plt.legend(fontsize=20.)    

    plt.show()
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
def simple_bend(F, a, L, EI, npts=100):
    """
        Computes the shape of a beam loaded with perpendicular force and simply supported (on both ends)

    Parameters
    ----------
    F : float
        The force applied to the beam.
    a : float
        Force application position from the first end of the beam.
    L : float
        The beam's length.
    EI : float
        The product of Young's Modulus & Quadratic moment of the beam's cross-setion.
    npts : int, optional
        Number of subdivions of calculation. (default 100)

    Returns
    -------
    x : numpy.ndarray
        First axis coordinates.
    y : numpy.ndarray
        Second axis coordinates.

    """
    ##########################################################
    # The differential equations that rules the deflection is:
    #           M.E.I = d²y/dx²
    ########################################################## 
    # Init Vars
    x = np.linspace(0,L,num=npts,dtype=np.float64)
    y = np.empty(npts)
    # Retrieve force position 'a'
    pos = 0
    for i in range(npts):
        if (x[i] > a):
            pos = i
            break
    # Dividing the work into two parts
    x1 = x[:pos]                    ; x2 = x[pos:]
    b = L - a ; L2 = L*L ; b2 = b*b ; x13 = x1*x1*x1 ; x23 = x2*x2*x2
    y1 = -x13 + (L*L - b*b)*x1      ; y2 = -x23 + L*(x2 - a)**3/b + (L*L - b*b)*x2
    # Joining the two parts
    fac = -F*b/(6*EI*L)
    y[:pos] = fac*y1
    y[pos:] = fac*y2
    
    return x, y
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