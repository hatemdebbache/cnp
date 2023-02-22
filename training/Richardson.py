import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from numpy.polynomial import Polynomial


def dh(f, h, x):
    '''
    Input:
        f: np.polynomial.Polynonimial type data. 
        h: floating point data.
        x: np.array type data.
    Output:
        return np.array type data of slope at each point x.
    '''
    a=f(x+h)
    b=f(x-h)
    dh=((a-b)/(2*h))
    return dh


def dh1(f, h, x):
    '''
    Input:
        f: np.polynomial.Polynonimial type data. 
        h: floating point data.
        x: np.array type data.
    Output:
        return np.array type data of slope at each point x.
    '''

    c=4*dh(f,h/2,x)
    d=dh(f,h,x)
    dh1=(c-d)/3
    return dh1


def error(f, hs, x_i):
    '''
    Input:
        f  : np.polynomial.Polynonimial type data. 
        hs : np.array type data. list of h.
        x_i: floating point data. single value of x.
    Output:
        return two np.array type data of errors by two methods..
    '''

    f_prime = f.deriv(1)
    Y_actual = f_prime(x_i)

    diff_error = []
    diff2_error = []

    for h in hs:
        # for each values of hs calculate the error using both methods
        # and append those values into diff_error and diff2_error list.

        # write your code here
      E1=Y_actual-dh(f,h,x_i)
      E2=Y_actual-dh1(f,h,x_i)
      diff_error.append(E1)
      diff2_error.append(E2)

      pass # delete this line
    
    print(pd.DataFrame({"h": hs, "Diff": diff_error, "Diff2": diff2_error}))

    return diff_error, diff2_error


#function to draw the actual function
def draw_graph(f, ax, domain=[-10, 10], label=None):
    data = f.linspace(domain=domain)
    ax.plot(data[0], data[1], label='Function')
    
    
fig, ax = plt.subplots()
ax.axhline(y=0, color='k')

p = Polynomial([2.0, 1.0, -6.0, -2.0, 2.5, 1.0])
p_prime = p.deriv(1)
draw_graph(p, ax, [-2.4, 1.5], 'Function')
draw_graph(p_prime, ax, [-2.4, 1.5], 'Derivative')

ax.legend()


fig, ax = plt.subplots()
ax.axhline(y=0, color='k')
hs = np.array([1., 0.55, 0.3, .17, 0.1, 0.055, 0.03, 0.017, 0.01])
e1, e2 = error(p, hs, 2.0)
ax.plot(hs, e1, label='e1')
ax.plot(hs, e2, label='e2')

ax.legend()