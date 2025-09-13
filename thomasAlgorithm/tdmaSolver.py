# tdmaSolver.py
# A Python implementation of the Thomas algorithm for solving tridiagonal systems of equations.
# Author: Jesse Blankenship
# Last Updated: 9/13/2025

import numpy as np

def tdmaSolver(a, b, c, d):
    """
    Solves the tridiagonal system of equations Ax = d using the Thomas algorithm.
    
    Parameters:
    a (np.array): the a coefficients of the main diagonal (phi at i)
    b (np.array): the b coefficients of phi at i+1
    c (np.array): the c coefficients of phi at i-1
    d (np.array): all costants for the equation defining phi at i
    """

    nPoints = np.size(a)
    
    # Initialize arrays
    P = np.zeros(nPoints)
    Q = np.zeros(nPoints)
    phi = np.zeros(nPoints)

    P[0] = b[0]/a[0]
    Q[0] = d[0]/a[0]
    for i in range(1, nPoints):
        P[i] = b[i]/(a[i]-c[i]*P[i-1])
        Q[i] = (d[i] + c[i]*Q[i-1])/(a[i] - c[i]*P[i-1])
    
    phi[-1] = Q[-1]
    for i in range(nPoints-2, -1, -1):
        phi[i] = P[i]*phi[i+1] + Q[i]

    return phi