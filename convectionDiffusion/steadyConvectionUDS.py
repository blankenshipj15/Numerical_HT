# steadyConvectionUDS is a solution implementaion for a 2D steady convection problem with known velocity and Dirichlet BCs.
# The problem is discretized using UDS first and then a deferred correction method using QUICK scheme to improve discretization order
# Author: Jesse Blankenship
# Last Updated: 10/16/2025

import numpy as np
import matplotlib.pyplot as plt

# Problem parameters
u = 1 # velocity in i direction
v = 1 # velocity in j direction
rho = 1 # density parameter
phiBottom = 1.0 # BC for bottom of domain
phiLeft = 0.0 # BC for left of domain
L = 1.0 # side length of square domain
nX = 80 # number of cells in i
nY = 80 # number of cells in j
deltaX = L / (nX-1)
deltaY = L / (nY-1)
maxIter = 10
tolerance = 1e-3

# Initialize solution matrix
phi = np.zeros((nY,nX))
# Solve for phi distribution using finite volume upwind difference scheme (first order)
aP = rho*u*deltaY + rho*v*deltaX # known coefficients from discretization
aW = rho*u*deltaY
aS = rho*v*deltaX
for i in range(nX):
    if i == 0: # first column of data (left edge)
        for j in range(nY):
            if j == 0:
                phi[j,i] = (aW*phiLeft + aS*phiBottom) / aP
            elif j < nY - 1: # all internal on left edge
                phi[j,i] = (aW*phiLeft + aS*phi[j-1,i]) / aP # upwinding to the cell below
            else:
                phi[j,i] = phi[j-1,i] # top cell is upwinded to interior value
    elif i < nX - 1:
        for j in range(nY):
            if j == 0:
                phi[j,i] = (aW*phi[j,i-1] + aS*phiBottom) / aP
            elif j < nY - 1: # all internal
                phi[j,i] = (aW*phi[j,i-1] + aS*phi[j-1,i]) / aP # upwinding to the cell below
            else:
                phi[j,i] = phi[j-1,i] # top cell is upwinded to interior value
    else:
        for j in range(nY):
            if j == 0:
                phi[j,i] = (aW*phi[j,i-1] + aS*phiBottom) / aP
            elif j < nY - 1: # all internal on right edge
                phi[j,i] = (aW*phi[j,i-1] + aS*phi[j-1,i]) / aP # upwinding to the cell below
            else:
                phi[j,i] = phi[j-1,i] # top cell is upwinded to interior value

# Function used for solving tridiagonal matrices (pulled from TDMA2D class)
def tdmaSolver(a, b, c, d):
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

# Iterate using deferred correction to recover a third order QUICK solution
rmsE = 1
a = np.zeros(nY); b = np.zeros(nY); c = np.zeros(nY); d = np.zeros(nY)
phiQuickE = np.zeros((nY,nX))
phiQuickW = np.zeros((nY,nX))
phiQuickN = np.zeros((nY,nX))
phiQuickS = np.zeros((nY,nX))
phiNew = np.zeros((nY,nX))
phiOld = phi.copy()

for iter in range(maxIter):
    # calculate the QUICK stencils for the interior cells
    # this needs updated every iteration using previous iteration phi values
    phiNew = np.zeros((nY,nX))
    phiQuickE = np.zeros((nY,nX))
    phiQuickW = np.zeros((nY,nX))
    phiQuickN = np.zeros((nY,nX))
    phiQuickS = np.zeros((nY,nX))
    for j in range(2, nY - 2):       
        for i in range(2, nX - 2):
            # using W,P,E
            phiQuickE[j,i] = phiOld[j,i] + (phiOld[j,i+1] - phiOld[j,i-1])/4.0 + (phiOld[j,i+1] + phiOld[j,i-1] - 2.0*phiOld[j,i])/8.0

            # using WW,W,P
            phiQuickW[j, i] = phiOld[j, i-1] + (phiOld[j, i] - phiOld[j, i-2]) / 4.0 + (phiOld[j, i] + phiOld[j, i-2] - 2.0 * phiOld[j, i-1]) / 8.0

            # using S,P,N
            phiQuickN[j, i] = phiOld[j, i] + (phiOld[j+1, i] - phiOld[j-1, i]) / 4.0 + (phiOld[j+1, i] + phiOld[j-1, i] - 2.0 * phiOld[j, i]) / 8.0

            # using SS,S,P
            phiQuickS[j, i] = phiOld[j-1, i] + (phiOld[j, i] - phiOld[j-2, i]) / 4.0 + (phiOld[j, i] + phiOld[j-2, i] - 2.0 * phiOld[j-1, i]) / 8.0

    for i in range(nX):
        if i == 0:
            for j in range(nY):
                a[j] = rho*u*deltaY + rho*v*deltaX
                b[j] = 0
                c[j] = 0 if j == 0 else rho*v*deltaX
                d[j] = rho*u*deltaY*phiLeft
                if j == 0:
                    d[j] += rho*v*deltaX*phiBottom
            phiNew[:, i] = tdmaSolver(a, b, c, d)

        elif i == nX - 1:
            for j in range(nY):
                a[j] = rho*u*deltaY + rho*v*deltaX
                b[j] = 0
                c[j] = 0 if j == 0 else rho*v*deltaX
                d[j] = rho*u*deltaY*phiOld[j, i-1]
                if j == 0:
                    d[j] += rho*v*deltaX*phiBottom
            phiNew[:, i] = tdmaSolver(a, b, c, d)

        else:  # interior vertical lines
            for j in range(nY):
                a[j] = rho*u*deltaY + rho*v*deltaX
                b[j] = 0
                c[j] = 0 if j == 0 else rho*v*deltaX
                d[j] = rho*u*deltaY*phiOld[j, i-1]
                if j == 0:
                    d[j] += rho*v*deltaX*phiBottom

                # Apply deferred correction
                if 2 <= i <= nX - 3 and 2 <= j <= nY - 3:
                    d[j] -= rho*u*deltaY * (phiQuickE[j, i] - phiOld[j, i]) # old value is upwinded from face e => P
                    d[j] += rho*u*deltaY * (phiQuickW[j, i] - phiOld[j, i-1]) # old value is upwinded from face w => W
                    d[j] -= rho*v*deltaX * (phiQuickN[j, i] - phiOld[j, i]) # old value is upwinded from face n => P
                    d[j] += rho*v*deltaX * (phiQuickS[j, i] - phiOld[j-1, i]) # old value is upwinded from face s => S

            phiNew[:, i] = tdmaSolver(a, b, c, d)

    # check for convergence
    print(iter)
    print(np.max(np.abs(phiNew-phiOld)))
    if np.max(np.abs(phiNew-phiOld)) <= tolerance:
        break
    phiOld = phiNew.copy()


plt.figure(1)
X = np.linspace(0, L, nX)
Y = np.linspace(0, L, nY)
X, Y = np.meshgrid(X, Y)
plt.contourf(X, Y, phiNew, levels=50, cmap='plasma')
plt.colorbar(label = "$\\phi$")
plt.title("$\\phi$ Distribution")
plt.xlabel("x/L")
plt.ylabel("y/L")
plt.show()

plt.figure(2)
yRangeAnalytic = np.linspace(0,L,1000)
analyticPhi = np.where(yRangeAnalytic < 0.5, 1, 0)
yRange = np.linspace(0,L,nY)
plt.plot(yRange,phi[:,nX//2])
plt.plot(yRange,phiNew[:,nX//2])
plt.plot(yRangeAnalytic,analyticPhi)
plt.legend(["UDS", "QUICK", "Analytic"])
plt.title("Vertical Centerline $\\phi$ Distribution")
plt.ylabel("$\\phi$")
plt.xlabel("y/L")
plt.show()

