# -*- coding: utf-8 -*-
"""
            Module: derive
            
            approximation numérique de dérivées
"""
import numpy as np

# ==============================================================
#                       Calcul de f'(x)
# ==============================================================
def dfp1(f, x, h):                                              # Ord = 1, DV = 1, C = 2
    return (f(x+h) - f(x))/h
# ==============================================================
def dfr1(f, x, h):                                              # Ord = 1, DV = 1, C = 2
    return (f(x) - f(x-h))/h
# ==============================================================
def dfc2(f, x, h):                                              # Ord = 2, DV = 2, C = 2
    """ 
        Estimation de la première dérivée de f en x 
    """
    return 0.5*(f(x+h) - f(x-h))/h
# ==============================================================
def Richardson1(f, x, hi, e, cmax, p):
    """
        Extrapolation de Richardson de la première dérivée de f en x <==> f'(x)
        
    Parameters
    ----------
    f : fonction à dériver
    x : point de dérivation
    hi : le pas initiale de calcule (sera divisé par deux à chaque itération)
    e : l'erreur sur l'estimation recherchée
    cmax : le coût maximale à prendre
    p (nb de colonne): limite d'ordre de pas 

    Returns
    -------
    c : le coût d'évaluation de la dernière valeur du tableau D[n-1,p-1]
    D : Le tableau de Richardson pour DFC2
    ee : l'erreur sur la dernière estimation
        
    N.B: La valeur la plus proche est D[n-1,p-1]
    That's the one you need!!!

    """
    n = int(0.5*(cmax - 3)) # On peut prédire le nombre de ligne selon la valeur de cmax
    
    # Initialisation des paramètres de retoure
    D = np.zeros((n, p))
    D[0, 0] = dfc2(f, x, hi)
    c = 3
    ee = float('inf')
    
    for i in range(1, n):
        hi *= .5
        D[i, 0] = dfc2(f, x, hi)
        c += 2
        for j in range(1, p):
            if (i >= j):
                A, B =  D[i,j-1], D[i-1, j-1]
                D[i, j] = A + (A - B)/(2**(2*j) - 1)
        
        ee = abs(D[i,p-1] - D[i-1,p-1]) if (i>=p) else ee
        if (ee < e):
            break

    return c, D[:i+1], ee
# ==============================================================
def dfp2(f, x, h):                                              # Ord = 2, DV = 2, C = 3
    return 0.5*(-3*f(x) + 4*f(x+h) - f(x+2*h))/h
# ==============================================================
def dfr2(f, x, h):                                              # Ord = 2, DV = 2, C = 3
    return 0.5*(f(x-2*h) -4*f(x-h) + 3*f(x))/h
# ==============================================================
#                       Calcul de f''(x)
# ==============================================================
def dfc22(f, x, h):                                             # Ord = 2, DV = 3, C = 3
    """ 
        Estimation de la dérivée deuxième de f en x 
    """    
    return (f(x+h) - 2*f(x) + f(x-h))/(h*h)
# =============================================================
def Richardson2(f, x, hi, e, cmax, p):
    """
        Extrapolation de Richardson de la dérivée seconde de f en x <==> f''(x)
        
    Parameters
    ----------
    f : fonction à dériver
    x : point de dérivation
    hi : le pas initiale de calcule (sera divisé par deux à chaque itération)
    e : l'erreur sur l'estimation recherchée
    cmax : le coût maximale à prendre
    p (nb de colonne): limite d'ordre de pas 

    Returns
    -------
    c : le coût d'évaluation de la dernière valeur du tableau D[n-1,p-1]
    D : Le tableau de Richardson pour DFC22
    ee : l'erreur sur la dernière estimation
        
    N.B: La valeur la plus proche est D[n-1,p-1]
    That's the one you need!!!

    """
    n = int(0.5*(cmax - 3)) # On peut prédire le nombre de ligne selon la valeur de cmax
    
    # Initialisation des paramètres de retoure
    D = np.zeros((n, p))
    D[0, 0] = dfc22(f, x, hi)
    c = 3
    ee = np.zeros(n)
    ee[0] = float('inf')
    
    for i in range(1, n):
        hi *= .5
        D[i, 0] = dfc22(f, x, hi)
        c += 2
        for j in range(1, p):
            if (i >= j):
                A, B =  D[i,j-1], D[i-1, j-1]
                D[i, j] = A + (A - B)/(2**(2*j) - 1)
        
        ee[i] = abs(D[i,p-1] - D[i-1,p-1]) if (i>=p) else ee[0]
        if (ee[i] < e):
            break

    return c, D[:i+1], ee
# ==============================================================
#                       Calcul de f'''(x)     
# ==============================================================
def dfc23(f, x, h):                                             # Ord = 2, DV = 4, C = 4
    """ 
        Estimation de la dérivée troisième de f en x 
    """
    return 2*(dfc2(f, x, 2*h) - dfc2(f, x, h))/(h*h)
# ==============================================================