import numpy as np


# =============================================================================
def comp_Trapeze(f, a, b, h):
    Nsub = (b - a)/h
    S = 0.5*(f(a) + f(b))
    for i in range(1, int(Nsub)):
        S = S + f(a + i*h)
    S = h*S
    
    return S
# =============================================================================
def Tnew(f, a, b, Nsub, Told):
        # Utilisé pour trouver la prochaine valeur de l'integral avec la methode des trapèzes
        # avec un pas double de Nsub, pour minimser le coût de calcule
        S = 0.
        H = (b - a)/Nsub
        for i in range(Nsub):
            S += f(a + (i + 0.5)*H)
    
        H = H/2.        
        
        return 0.5*Told + H*S
# =============================================================================
def dfc2(f, x, h):                                              # Ord = 2, DV = 2, C = 2
    return 0.5*(f(x+h) - f(x-h))/h
# =============================================================================
def Richardson(func, pas, e, n, p, **kwargs):
    """
        Apply the Richardson extrapolation for a numerical method represented by 'func'

    Parameters
    ----------
    func : the procedure for the numerical method
    pas : the initial step of calculation
    n (rows): number of steps (each step is divided by two) 
    p (columns): number of orders (each order is multiplied by two) 
    *args : arguments for the numerical method, passed 'as is' to func

    Returns
    -------
    c : cost of calculation (number of evaluations of the function given to the numerical method)
    Ret : Richardson's extrapolation table
    ee : error estimation of returned value

    """
    
    # Initialising returned variables
    c = 0
    Ret = np.zeros((n, p))
    ee = float('inf')
    
    Ret[0,0] = func(h=pas, **kwargs)
    # f, x = kwargs['f'], kwargs['x']
    # Ret[0,0] = 0.5*(f(x+pas) - f(x-pas))/pas
    i = 0
    while(ee > e and i < n):
        pas *= .5
        Ret[i, 0] = func(h=pas, **kwargs)
        c += 1
        for j in range(1, p):
            if (i >= j):
                A, B =  Ret[i,j-1], Ret[i-1, j-1]
                Ret[i, j] = A + (A - B)/(2**(2*j) - 1)
        ee = abs(Ret[i,p-1] - Ret[i-1,p-1]) if (i > p) else float('inf')
        i += 1
        
    return c, Ret[:i], ee

from math import sqrt

def f(x):
    return sqrt(x)

c, ret, ee = Richardson(comp_Trapeze, pas=0.1024, e=1e-15, n=10, p=3, f=f, a = 2, b = 5)
print(ret)