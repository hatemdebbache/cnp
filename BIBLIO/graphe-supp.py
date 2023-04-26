"""
    Module : Graphe Supplémentaire
    
    Fonctions supplémentaires pour traçage des fonctions/nuages de points.
    
    By: HD 
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler
import numpy as np
import edo

mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
plt.style.use('dark_background')
# plt.style.use('classic')

# ==============================================================
def TracePortrait(f, a, b, PtCs, npts, h, r=0.1, ntraj=20, methode = 'RK', relie=True, couleur = 'red', epaisseur = 1.0):
    """
        Tracer le portrait de phase d'un système d'équations différentielles ordinaires autonomes d'ordre 2.
    
        Cette fonction prend des points des cercles centrés en chaque point critique *PtCs* et de rayon *r* comme condition initial d'intégration.
    
    Parameters
    ----------
    f : fonction
        A function that defines the system of ODEs as (y', z') = f(t, y, z).
    a : float
        The left endpoint of the interval.
    b : float
        The right endpoint of the interval.
    PtCs : ndarray
        List of critical points coordinates (ex: [(0,0), (1,1)...] ).
    npts : int
        Number of points to be plotted per trajectory.
    h : float
        Integration step.
    r : float, optional
        The radius of CI cercle. The default is 0.1.
    ntraj : int, optional
        Number of trajectories per critical point. The default is 20.
    methode : string, optional
        The numerical solver to be used. The default is 'RK'.
    relie : bool, optional
        Wether or not to link points (plot/scatter). The default is True.
    couleur : string, optional
        Color of lines. The default is red.
    epaisseur : float, optional
        Linewidth of plot. The default is 1.0px.
        
    List of available numerical solvers
    ----------
    *Euler* : Euler's explicit method (see edo.Euler2).
    
    *Heun* : modified Euler's method (see edo.Heun2).
    
    *PM* : Midpoint method (see edo.PointMilieu2).
    
    *RK* : Runge-Kutta fourth order method (see edo.RungeKutta2).

    """
    Solvers = {
        'Euler' : edo.Euler2,
        'Heun' : edo.Heun2,
        'PM' : edo.PointMilieu2,
        'RK':edo.RungeKutta2,
        }
    solver = Solvers[methode]
    
    m = int((b - a)/h + 0.5)
    theta = np.linspace(0,2*np.pi,num=ntraj,endpoint=False)
    hor, ver = r*np.cos(theta), r*np.sin(theta)
    for pc in PtCs:
        Ya = pc[0] + hor ; Za = pc[1] + ver
        for yo, zo in zip(Ya, Za):
            c,t,y,z = solver(f, a, b, yo, zo, npts, m)
            if (relie):
                plt.plot(y,z, color = couleur, lw = epaisseur)
            else:
                plt.scatter(y,z, color = couleur, lw = epaisseur)
# ==============================================================
def TraceFlux(f, ymin, ymax, zmin, zmax, step):
    """
        Draw the flow field of autonomous ODE system.

    Parameters
    ----------
    f : TYPE
        DESCRIPTION.
    ymin : TYPE
        DESCRIPTION.
    ymax : TYPE
        DESCRIPTION.
    zmin : TYPE
        DESCRIPTION.
    zmax : TYPE
        DESCRIPTION.
    step : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """