# tdma1Dtests is the main file used for testing the TDMA function in a 1D heat conduction problem.
# Author: Jesse Blankenship
# Last Updated: 9/13/2025

import numpy as np
import matplotlib.pyplot as plt
from tdmaSolver import tdmaSolver

# Problem parameters
gamma = 400.0               # conductivity
tBcLeft1 = 500.0            # Dirichlet left BC (Case 1 & 2)
tBcRight1 = 300.0           # Dirichlet right BC (Case 1)
fluxBcRight2 = 3000.0       # Neumann right BC (Case 2)
domainLen = 6.0             # total length of 1D domain
nPoints = 3                 # number of cell centroids
deltaX = 2.0                # distance between centroids

# Case 1: Dirichlet BCs
a1 = [gamma/deltaX + 2*gamma/deltaX,
      gamma/deltaX + gamma/deltaX,
      2*gamma/deltaX + gamma/deltaX]

b1 = [gamma/deltaX, gamma/deltaX, 0]
c1 = [0, gamma/deltaX, gamma/deltaX]
d1 = [2*gamma*tBcLeft1/deltaX, 0, 2*gamma*tBcRight1/deltaX]

T1 = tdmaSolver(a1, b1, c1, d1)

print("Case 1:")
print(tBcLeft1)
for i in range(nPoints):
    print(T1[i])
print(tBcRight1, "\n")

# Case 2: Dirichlet left, Neumann right
a2 = [gamma/deltaX + 2*gamma/deltaX,
      gamma/deltaX + gamma/deltaX,
      gamma/deltaX]

b2 = [gamma/deltaX, gamma/deltaX, 0]
c2 = [0, gamma/deltaX, gamma/deltaX]
d2 = [2*gamma*tBcLeft1/deltaX, 0, -fluxBcRight2]

T2 = tdmaSolver(a2, b2, c2, d2)

print("Case 2:")
print(tBcLeft1)
for i in range(nPoints):
    print(T2[i])
    if i == nPoints - 1:
        # Reconstruct wall temperature from Neumann BC
        ghost_T = (((2*gamma/deltaX)+(gamma/deltaX))*T2[-1] - (gamma/deltaX)*T2[-2]) / (2*gamma/deltaX)
        print(ghost_T, "\n")
