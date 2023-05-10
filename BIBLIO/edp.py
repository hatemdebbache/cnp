# -*- coding: utf-8 -*-
"""
        Module: EDP
        Resolution des equations differentielles aux derivees partielles.
        
"""
import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
def get_kernel(eq, dx, dy):
    """
    Generate the calculation kernel for Gauss-Seidel method to solve 2-D Partial Differential Equation.
    
    Parameters
    ----------
    eq : ndarray
        List of coefficients of differential operators that define the equation \
            (See Notes for more details).
    dx : float
        Step in x-axis.
    dy : float
        Step in y-axis.

    Returns
    -------
    kernel : ndarray
        Matrix representing the calculation kernel.
        
    Notes
    -----
    - The general form of a differential equation is:
        $$c_0 U_x + c_1 U_y + c_2 U_{xx} + c_3 U_{yy} + c_4 U_{xy} + c_5 U_{xxx} + ... = f(x)$$
    *Provide those coefficient in list and pass it as `eq`*.
    
    ``*Warning!*``
    For now it works for differential operators up-to Uxy and f(x) = 0, further updates will include the general case code.
    
    - This function uses the discretised form of the differential operators to build the final calculation kernel,\
        when using small integrating steps truncation errors may build up and significant error shows in the result.\
        Use with caution.
        
    - To solve the differential equation in a square grid of (n,m) elements, apply the formula:
        
    >>> for i in range(n):
    >>>     for j in range(m):
    >>>         U[i,j] = sum(sum(U[i-k:i+k+1, j-k:j+k+1]*kernel)) 
    
    *(The product of the matrix and kernel is element-wise).*
    
    `k` is half the kernel order which is dependent on the order of the differential equation.
    Ex: the kernel size of a second-order PDE is 3, k = 1.
    
    """
    Ux = 0.5/dx*np.array([[0., -1., 0.], 
                          [0., 0., 0.],
                          [0., 1., 0.]])
    Uy = 0.5/dy*np.array([[0., 0., 0.],
                          [-1., 0., 1.],
                          [0., 0., 0.]])
    Uxx = 1/(dx*dx)*np.array([[0., 1., 0.],
                              [0., 0., 0.],
                              [0., 1., 0.]]) 
    Uyy = 1/(dy*dy)*np.array([[0., 0., 0.],
                              [1., 0., 1.],
                              [0., 0., 0.]]) 
    Uxy = 0.25/(dx*dy)*np.array([[1., 0., -1.],
                                 [0., 0., 0.],
                                 [-1., 0., 1.]])
    
    kernels = [Ux, Uy, Uxx, Uyy, Uxy]
    cst = np.array([0., 0., 2/(dx*dx), 2/(dy*dy), 0.])
    
    kernel = np.zeros((3,3))
    for i in range(len(kernels)):
        kernel += eq[i]*kernels[i]
    
    kernel /= sum(eq*cst)
    
    return kernel
# ============================================================================
def heatmap(U):
    """
    Show the heatmap of an array

    Parameters
    ----------
    U : ndarray
        The array to be represented.
    """
    x = np.arange(U.shape[1]) ; y = np.arange(U.shape[0])
    
    plt.clf()
    plt.xticks(x) ; plt.yticks(y)

    plt.rcParams.update({'axes.grid' : True,
                         'grid.color' : 'white', 
                         'grid.linewidth' : 2.,
                         
                         'xtick.labeltop' : True,   
                         'xtick.labelbottom' : False,
                          
                         'ytick.labelleft' : True,   
                         'ytick.labelright' : False,   
                         })
    
    plt.imshow(U,cmap='hot')
    plt.colorbar()
    
    plt.show()