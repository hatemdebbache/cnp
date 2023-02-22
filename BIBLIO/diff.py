"""
        Module : Diff
        Méthodes pour résoudre des équations différentielles
"""
import numpy as np
import syslin

# ============================================================================
#
#           Problèms aux limites: cas linéaire
#
# ============================================================================
def EqDiffLin2(P,Q,R,a,b,Va,Vb,n): 
    """
        Résolution d'une équation différentielle de deuxième ordre de la forme:    
            y'' = P(x).y' + Q(x).y + R(x)
        Avec des conditions aux limites y(a) = Va, y(b) = Vb
        n : nombre de pas d'integration
    """
    # Initialisation    
    L = np.zeros(n)
    D = np.zeros(n)
    U = np.zeros(n)
    b = np.zeros(n) 
    # Pas d'integration
    h = (b - a)/n
    # Conditions aux limites
    b[0] = R(a) - (1 + 0.5*h*P(a))*Va
    b[n-1] = R(b) - (1 - 0.5*h*P(b))*Vb
    # Construction du système linéaire équiv
    for k in range(1,n-1):
        L[k] = 1 - 0.5*h*P(a + (k+1)*h)
        D[k] = -h*h*Q(a + k*h)
        U[k] = 1 + 0.5*h*P(a + (k-1)*h)
        b[k] = h*h*R(a + k*h)
    # Résolution
    y = syslin.Gauss_Thomas2(L,D,U,b)
    
    return y
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================
# ============================================================================