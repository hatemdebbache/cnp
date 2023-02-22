"""
    Matrices Module
    
    functions to deal with matrices and vectors for solving linear problems

    MADE BY: Hatem D.
    
    Note: these methods are functional but need improvement 
        (to avoid errors in critical conditions)
"""
#==================================================
import numpy as np
from numpy import linalg as LA
from math import sqrt, cos, sin, acos, atan2
import sys
sys.path.insert(0, "./../BIBLIO/")
import cbase
#==================================================

# -------
eps = np.finfo(np.float64).eps
#==================================================
def is_perp(u, v):
    """
    Parameters
    ----------
    u, v : 1-D arrays
    
    Returns
    -------
    Boolean: True if u and v are perpendicular, else False.
    """
    u, v = check_ndarrays(u), check_ndarrays(v)
    # Verify if vectors
    if u.ndim == 1 and v.ndim == 1:
        return not(u@v)
    else:
        print("Parameters should be 1-D array")
#==================================================
def module(u):
    """
    Calculates the length of vector u
    
    Parameters
    ----------
    u : 1-D array
    
    Returns
    -------
    float: 2-norm of NumPy.LinearAlgebra
    """
    u = check_ndarrays(u)
    return LA.norm(u)
#==================================================
def op_vec(u):
    """
    Operateur produit vectoriel

    Parameters
    ----------
    u : 1-D array
    
    Returns
    -------
    2-D array: the matrice for vectoriel product u ^ v = [*u].v
    """
    a  = np.zeros((3,3))
    a[0][1] = -u[2]
    a[0][2] = u[1]
    a[1][2] = -u[0]
    a = a - a.transpose()
    return check_ndarrays(a)
#==================================================
def check_ndarrays(*vcts):
    """
    Check if array are of type np.ndarray, else transform to corresponding type

    Parameters
    ----------
    *vcts : one or more list(/ndarray)
    
    Returns
    -------
    List of ndarray of the given parameters
    """
    result = []
    for v in vcts:
        if type(v) != np.ndarray:
            v = np.array(v)
        result.append(v)

    return tuple(result) if len(result) > 1 else result[0]
#==================================================
def identity(n = 3):
    """
    Get I the identity matrix of nth-order
    I.M = M.I = M 
    M is a matrix of nth-order

    Parameters
    ----------
    n : int, optional
        the size of the identity matrice. The default is 3.

    Returns
    -------
    I : np.ndarray
        Zeros matrice with ones in diagonal.
    """
    a = np.zeros((n,n))
    for i in range(n):
        a[i,i] = 1
    
    return a
#==================================================
def normalize(u):
    """
    Get the unit vector of same direction

    Parameters
    ----------
    u : 1-D array
        
    """
    u = check_ndarrays(u)
    return u / module(u)
#==================================================
def project(c, p1, p2):
    """
    Project the point C on the line (P1P2)

    Parameters
    ----------
    c : 1-D array
        Position vector of point C.
    p1 : 1-D array 
        Position of first point defining the line.
    p2 : 1-D array
        Position of second point defining the line.

    Returns
    -------
    h : 1-D array
        Position of H, projection of C on (P1P2).

    """
    c, p1, p2 = check_ndarrays(c, p1, p2)
    u = normalize(p2-p1)
    u2 = op_vec(u)
    h = (identity() + u2@u2)@c - (u2@u2)@p1
    
    return h
#==================================================
def cercline(c, r, p1, p2):
    """
    Check the intersection between the circle (C, R) and the line segment [P1P2]

    Parameters 
    ----------
    c : 1-D array
        Position vector of point C, center of the circle.
    r : float
        Radius of the cercle.
    p1 : 1-D array 
        Position vector of point P, first point defining the line segment.
    p2 : 1-D array
        Position vector of point P, second point defining the line segment.

    Returns
    -------
    Integer
        0 if there's no intersection between the cercle and the line segment.
        1 if there's intersection in one point (tangent or not).
        2 if there's intersection in two different points.

    """
    c, p1, p2 = check_ndarrays(c, p1, p2)
    d = int(module(c-p1) > r)
    d += int(module(c-p2) > r)
    if d > 1:
        h = project(c, p1, p2)
        d = module(h-c)
        if d < r:
            return 2
        elif d == r:
            return 1
        else:
            return 0
    else:
        return d
# ==================================================
def circcirc(x1, y1, r1, x2, y2, r2):
    r3 = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    indx1 = r3 > r1 + r2
    indx2 = r2 > r3 + r1
    indx3 = r1 > r3 + r2
    indx4 = r3 < 10*eps and (r2 - r1) < 10*eps
    indx = indx1 + indx2 + indx3 + indx4
    if indx > 0:
        return None
    
    anought = atan2((y2 - y1), (x2 - x1))
    
    aone = acos(-(r2*r2 - r1*r1 - r3*r3)/(2*r1*r3))
    
    alpha1 = anought + aone
    alpha2 = anought - aone
    
    p1 = [x1 + r1*cos(alpha1), y1 + r1*sin(alpha1)]
    p2 = [x1 + r1*cos(alpha2), y1 + r1*sin(alpha2)]
    
    return p1, p2
# ==================================================
def linecirc2(x1, y1, x2, y2, xc, yc, r):
    slope = x2 - x1
    if (slope < 10*eps):
        # Infinite slope case
        if (abs(xc - x1) > r):
            x, y = None, None
        else:
            step = sqrt(r*r - (xc - x1)*(xc - x1))
            x = [x1, yc + step]
            y = [x1, yc - step]
    
    else:
        slope = (y2 - y1)/slope
        intercpt = y1 - slope*x1
        
        a = 1 + slope*slope
        b = 2*(slope*(intercpt - yc) - xc)
        c = yc*yc + xc*xc + intercpt*intercpt - 2*yc*intercpt - r*r
        
        s, w1, w2 = cbase.esd(b/a, c/a)
        if s == -1:
            x, y = None, None
        else:
            x = [w1, intercpt + slope*w1]
            y = [w2, intercpt + slope*w2]
    
    
    return x, y
# ==================================================
def linecirc(slope, intercpt, xc, yc, r):
    if (abs(slope) < 10*eps):
        # Infinite slope case
        if (abs(xc - intercpt) > r):
            x, y = None, None
        else:
            step = sqrt(r*r - (xc - intercpt)*(xc - intercpt))
            x = [intercpt, yc + step]
            y = [intercpt, yc - step]
    
    else:
        a = 1 + slope*slope
        b = 2*(slope*(intercpt - yc) - xc)
        c = yc*yc + xc*xc + intercpt*intercpt - 2*yc*intercpt - r*r
        
        s, w1, w2 = cbase.esd(b/a, c/a)
        if s == -1:
            x, y = None, None
        else:
            x = [w1, intercpt + slope*w1]
            y = [w2, intercpt + slope*w2]

    return x, y