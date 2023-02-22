"""
    Module: cbase
    
    Des fonctions élémentaires
"""
import numpy as np
from math import sqrt, tanh
# ==============================================================
def signe(x):
    """
        Return
        -------
        0 : x est nulle
        1 : x est positif
        -1 : x est négatif
    """
    if(x < 0):
        y = -1
    elif(x > 0):
        y = 1
    else:
        y = 0
    
    return y
# ==============================================================
def signe2(x):
    """
        Return
        -------
        0 : x est nulle
        1 : x est positif
        -1 : x est négatif
    """
    return tanh(1e9*x)
# ==============================================================
def Heaviside(x):
    return 0.5*(1+signe2(x))
# ==============================================================
def intersection_cercle(c, r, p1, p2):
    """
        Check the intersection status of a circle of center c (1-D array with coordinates) and rayon r
            with a line segment between the points p1 and p2 (1-D arrays with coordinates)
    
        Return:
            0   The segment doesn't intersect with the circle (or the segment intersects with the circle in one point, 
                    the other point is out of domain)
            1   The line is tangent to the circle
            2   The segment intersects with the circle in two different points
    """
    xc, yc = c
    x1, y1 = p1
    x2, y2 = p2
    
    a = (x2-x1) 
    b = (y2-y1) 
    c = (xc+x1) 
    d = (yc+y1)
    
    a1 = 2*(a*c + b*d)/(a**2 + b**2)
    a2 = (c**2 + d**2 - r**2)/(a**2 + b**2)
    
    s, x, y = esd(a1, a2)
    if s == 1:
        if (x >= 0 and x <= 1) and (y >= 0 and y <= 1):
            return 2
        return 0
    else:
        return 1 + s
# ==============================================================
def ESD(b,c):
    # Calcule du determinant
    det = b*b - 4*c
    # Ctrl des cas possibles
    s = signe(det)
    if (s == 1):
        root = sqrt(det)
        x = 0.5*(-b + root)
        y = 0.5*(-b - root)
        if(abs(x) > abs(y)):
            x,y = y,x
    elif (s == 0):
        x = -0.5*b
        y = 0
    else:
        x = -0.5*b
        y = 0.5*sqrt(abs(det))
    
    return s, x, y
# ==============================================================
def esd(b,c):
    """ 
        Résolution d'équation de seconde degrée de la forme x*x + b*x + c = 0
        (amélioré)
        
        Return
        -------
        s : signe(déterminant)
        Si s >= 0:
            x, y : racines de l'équation
        Sinon:
            x : partie réel
            y : partie imaginaire
    """
    x = -0.5*b
    y = x*x - c
    s = signe(y)
    y = sqrt(s*y)
    if (s == 1):
        if (x > 0):
            y = x + y
        else:
            y = x - y
            
        x = c/y
    
    return s,x,y
# ==============================================================
def Poly_Lagrange(X, Y, x):
    """
        Interpolation avec polynôme de Lagrange
        
        Return
        -------
        y = P(x) : l'image de x par le polynôme d'interpolation des points (Xi, Yi)
    """
    N = len(X)
    y = 0.
    for i in range(0, N):
        p = 1.
        for j in range(0, N):
            if not (i == j):
                p = p*(x - X[j])/(X[i] - X[j])
        y = y + p*Y[i]
    
    return y
# ==============================================================
def Poly_Legendre(n, x):
    if (n == 0):
        return 1
    if (n == 1):
        return x
    
    P0 = 1
    P1 = x
    for i in range(1,n): # Verify the formula 
        P = ((2*i + 1)*x*P1 - i*P0)/(1+i)
        P0, P1 = P1, P
    
    return P
# ==============================================================
def Somme_Legendre(n, C, x):
    if (n <= 1):
        return C[0]
    P0, P1 = 1, x
    S = C[0] + C[1]*P1
    for i in range(2,n): # Verify the formula 
        P = ((2*i - 1)*x*P1 - (i-1)*P0)/i
        S += C[i]*P
        P0, P1 = P1, P
    
    return S
# ==============================================================
def Horner(C, x):
    """
        Evaluation éfficace de P(x) 'Méthode de Horner'
        
        Paramètres
        -----------
        C : list des coefficient du polynôme C = [Co, C1, C2 ... Cn] => P(x) = Cn*x^n + Cn-1*x^(n-1)... + C1*x + Co
        x : point ou on veut evaluer P(x)
        
        Return
        -----------
        y = P(x)
    """
    N = len(C)
    i = N - 1
    y = C[i]
    while (i > 0):
        i = i - 1
        y = x*y + C[i]
    
    return y
# ==============================================================
def affiche_tableau(T, n, titre=""):
    """
        Afficher un tableau des réels avec n chiffre significatif de façon laisible

    Parameters
    ----------
    T : Le tableau des valeurs à afficher
    n : le nombre de chiffres significatif
    titre : on peut ajouter un titre au tableau
    
    """
    s = '###'
    s *= n
    print(f'{s}{titre}{s}')
    for line in T:
        for cell in line:
            s = f'  {cell:.{n}g}'
            while(len(s) < n+10):
                s += ' '
            print(s, end='#')
        print()
    s = '###'
    s *= n
    print(f'{s}{s}')
# ==============================================================
def check_ndarray(*args, silent=True):
    """
    Check if given arguments are of type 'numpy.ndarray', else convert them.

    Parameters
    ----------
    *args : any_type
        The objects to be verified.
    silent : boolean, optional
        Do not to display the message 'converting array'.
        
    Raises
    ------
    TypeError
        One of the obejcts are not array_like.
    
    Returns
    -------
    A : numpy.ndarray
        The ndarray object.
    """
    ret = []
    for A in args:
        if not isinstance(A, np.ndarray):
            try:
                if not silent: print(f"Converting argument {args.index(A)} to ndarray...")
                A = np.array(A, dtype=np.float64)
            except:
                raise TypeError(f"Argument in position {args.index(A)} cannot be converted to ndarray")
        ret.append(A)
    return ret[0] if (len(ret) <= 1) else tuple(ret)
# ==============================================================