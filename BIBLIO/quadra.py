# -*- coding: utf-8 -*-
"""
                Module Quadra
            Approximation des intégrales
"""
import numpy as np

# ===================================================
#                   Formules Composites
# ===================================================
def comp_PointMilieu(f, a, b, Nsub):
    H = (b - a)/Nsub
    S = 0.
    for i in range(0, Nsub):
        S = S + f(a + (i + 0.5)*H)
    S = H*S
    
    return S
# ===================================================
def comp_Trapeze(f, a, b, Nsub):
    H = (b - a)/Nsub
    S = 0.5*(f(a) + f(b))
    for i in range(1, Nsub):
        S = S + f(a + i*H)
    S = H*S
    
    return S
# ===================================================
def comp_Simpson(f, a, b, Nsub):
    S1 = comp_Trapeze(f, a, b, Nsub)
    S2 = comp_PointMilieu(f, a, b, Nsub)
    S = (S1 + 2.*S2)/3.
    
    return S
# ===================================================
def comp_Simpson2(f, a, b, Nsub):
    T = comp_Trapeze(f, a, b, Nsub)
    H = (b - a)/Nsub
    H6 = H/6
    S = 0.
    for i in range(0, Nsub):
        c = a + (i + 0.5)*H
        S = S + f(c + H6) + f(c - H6)
    S = 0.25*(1.5*H*S + T)
    
    return S
# ===================================================
#                   Méthodes à pas variable
# ===================================================
def PV_trapeze(f, a, b, e, cmax):
    """
        Estimation de l'integrale de f dans [a, b] avec la méthode des Trapeze à pas variable
    
    Parameters
    ----------
    f : fonction à intégrer
    
    a, b : domaine d'integration
    
    e : l'erreur sur l'approximation recherchée
    
    cmax : le coût maximale de calcule permis
    
    Returns
    -------
    c : le coût consomé par la méthode (0 si cmax est atteint)
    
    T : la plus prôche estimation de la valeur réelle
    
    ee : l'estimation de l'erreur sur T
    """
    h = 0.5*(b - a)
    Tanc = h*(f(a) + f(b))
    T = 0.5*Tanc + h*f(a + h)
    ee = abs(T - Tanc)
    c = 3
    m = 1
    while (True):
        m = 2*m
        c = c + m
        if (c > cmax):
            c = 0
            break
        S = 0.        
        for i in range(0, m): 
            S = S + f(a + (i + 0.5)*h)
        h = 0.5*h
        Tanc = T
        T = 0.5*Tanc + h*S
        ee = abs(T - Tanc)
        if (ee < e):
            break
    
    return c, T, ee
# ===================================================
def Romberg(f, a, b, Nsub, e, n, p):
    """
        Extrapolation de Richardson de l'intégrale de f dans [a, b] à l'aide de la méthode des Trapèzes
        
    Parameters
    ----------
    f : fonction à intégrer
    [a, b] : domaine de l'intégration
    Nsub : le nombre de subdivisions initiale de calcule (sera multiplié par deux à chaque itération)
    e : l'erreur sur l'estimation recherchée
    n (nb de ligne) : nombre de Nsub à prendre
    p (nb de colonne): limite d'ordre de pas 

    Returns
    -------
    c : le coût d'évaluation de la dernière valeur du tableau D[n-1,p-1]
    I : Le tableau de Richardson pour comp_Trapeze
    ee : l'erreur sur la dernière estimation
        
    N.B: La valeur la plus proche est D[n-1,p-1]
    That's the one you need!!!

    """
    # Initialisation des paramètres de retoure
    I = np.zeros((n, p))
    I[0, 0] = comp_Trapeze(f, a, b, Nsub)
    c = 2
    ee = np.zeros(n)
    ee[0] = float('inf')
    
    for i in range(1, n):
        c += Nsub
        Nsub *= 2
        I[i, 0] = comp_Trapeze(f, a, b, Nsub)
        for j in range(1, p):
            if (i >= j):
                A, B =  I[i,j-1], I[i-1, j-1]
                I[i, j] = A + (A - B)/(2**(2*j) - 1)
        
        ee[i] = abs(I[i,p-1] - I[i-1,p-1]) if (i>=p) else float('inf')
        if (ee[i] < e):
            break

    return c, I[:i+1], ee
# ===================================================