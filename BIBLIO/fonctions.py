# -*- coding: utf-8 -*-
"""
                            Module Fonctions
                Des procedures essentielles pour les fonction
"""
import sys 
sys.path.insert(0, 'E:/EMP/CNP/BIBLIO/')
import numpy as np
from math import cos, sin
from MatVec import vec_interval

# ==============================================================
def trouver_abscisse(f, a, pas = 0.1):
    """
        Return the first (approx.) point x where f(x) = 0
        
        approxi.: only equidistant points are considered with intervals of error of 'pas'
    """
    x = a
    n = 2
    a = f(x)
    b = f(x+pas)
    while(a*b >= 0):
        n += 1
        x += pas
        a = b
        b = f(x+pas)

    return x + 0.5*pas, n
# ==============================================================
def nb_nul(f, a, b, pas = 0.1):
    """ 
        Return the (approx.) number of occurences of f(x) = 0
        
        approxi.: only equidistant points are considered with intervals of error of 'pas'
    """
    x = a
    n = 0
    while(x < b):
        x = trouver_abscisse(f, x)
        n = n + 1
    
    return n
# ==============================================================
def create_circle(theta, r, xo, yo):
    """ 
        Return the coordinates (x,y) of the arc of center (xo, yo), rayon r, and angle theta
    """
    n = len(theta)
    x = np.zeros(n)
    y = np.zeros(n)
    for i in range(n):
        x[i] = xo + r*cos(theta[i])
        y[i] = yo + r*sin(theta[i])
    return x, y
# ==============================================================
def discretise_fonction(f, xmin, xmax, npts):
    """ 
        Return list of function values for n point in interval [xmin, xmax]
    """
    x = vec_interval(xmin, xmax, npts)
    y = np.zeros(npts)
    for i in range(npts):
        y[i] = f(x[i])
    
    return x, y
# ==============================================================
def derive(g, n, x, npts, h):
    # BENEFIT: Estimate all derivatives (up to 'npts' derivative, lower orders are more precise) at once 
    # UNSTABLE FOR A RANGE OF VALUES FOR NPTS, H
    # Maximum npts of 113 for h 1e-3 before breaking
    
    # With optimum it works so much better
    
    # This method is still based on finite differences, look for alternative methods
    

    A = np.zeros((npts,npts))  # A holds coeffs of f derivatives in x, each line evaluated in a pt around x 
    
    c = int((npts - 1)/2)
    A[c, 0] = 1
    
    for i in range(c):
        f = 1
        for j in range(npts):
            A[i, j] = (h*(c - i))**j/f
            A[-i - 1, j] = (h*(i - c))**j/f
            f = f*(j+1)
    
    Y = np.zeros(npts)  # Y holds g(x+ih) for all pts used to make A
    Y[c] = g(x)
    for i in range(c):
        Y[i] = g(x + (c - i)*h)
        Y[-i - 1] = g(x + (i - c)*h)
    
    v = np.linalg.solve(A, Y) # v holds the value of all f derivatives in x

    return v[n]
# ==============================================================
def optimum(g, n, x):
    # this method find the best nb of pts used to estimate the n derivative of g in x
    
    # ADD METHOD FOR FINDING BEST (nb, h) 'as coupled values'
    nb = n + 3 if n%2 == 0 else n + 2
    
    def relative_error(g, n, x, nb):    
        w = derive(g, n, x, nb, .1)
        nb += 2
        w2 = derive(g, n, x, nb, .1)
        e = abs((w2 - w)/w)
    
        return e, w2
    
    e, w = relative_error(g, n, x, nb)
    nb += 2
    e2, w = relative_error(g, n, x, nb)
    
    while (e2 < e):
        e = e2
        nb += 2
        e2, w = relative_error(g, n, x, nb)
    
    return nb, w
# ==============================================================