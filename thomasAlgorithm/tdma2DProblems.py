# tdma2DProblems is a solution implementaion for a 2D steady state heat conduction problem with Dirichlet BCs and TDMA solver.
# Author: Jesse Blankenship
# Last Updated: 9/21/2025

import numpy as np
import matplotlib.pyplot as plt
from tdmaSolver import tdmaSolver
from TDMA2D import TDMA2D
import time

# Problem parameters
gamma = 5.0 # conductivity
tBcLeft = 20.0 # Dirichlet left BC
tBcRight = 100.0  # Dirichlet right BC
tBcTop = 100.0 # Dirichlet top BC
tBcBottom = 20.0 # Dirichlet bottom BC
height = 1.0  # domain height
width = 1.0  # domain width
nX = 50 # number of cell centroids in x
nY = 50 # number of cell centroids in y
deltaX = width / (nX - 1) # distance between centroids in x
deltaY = height / (nY - 1) # distance between centroids in y
maxIter = 1000
tolerance = 1e-4

# Initialize temperature array and coefficient arrays
T = np.zeros((nY, nX))
T_prev = np.mean(tBcTop+tBcBottom+tBcRight+tBcLeft) * np.ones((nY, nX))
# array sized with nY but this could also be nX (assuming square domain)
a = np.zeros((nY,1)); b = np.zeros((nY,1)); c = np.zeros((nY,1)); d = np.zeros((nY,1))

time_start = time.time()
# Initialize TDMA solver
T = TDMA2D(nX, nY, width, height, gamma, tolerance, maxIter, tBcLeft, tBcRight, tBcTop, tBcBottom)
# Solve using the solve method
T = T.solve()
time_end = time.time()
print("Time to converge: ", time_end - time_start)

# append the BCs the solution for plotting
T_numerical = np.zeros((nY+2, nX+2))
T_numerical[0,:] = tBcBottom
T_numerical[-1,:] = tBcTop
T_numerical[:,0] = tBcLeft
T_numerical[:,-1] = tBcRight
T_numerical[1:-1,1:-1] = T

TSliceHorizontal = T_numerical[nY//2,:]
TSliceVertical = T_numerical[:,nX//2]

# Analytic Solution for comparison
numModes = 200
T_A = np.zeros((nY+2, nX+2))
c_A = np.zeros(numModes)
lambda_A = np.zeros(numModes)
T_B = np.zeros((nY+2, nX+2))
theta = np.zeros((nY+2, nX+2))
T_analytic = np.zeros((nY+2, nX+2))
c_B = np.zeros(numModes)
lambda_B = np.zeros(numModes)
x = np.linspace(0, width, nX+2)
y = np.linspace(0, height, nY+2)

# solution is referenced from Heat Transfer, Nellis and Klein section 2.4.2
for i in range(numModes):
    lambda_A[i] = (i+1) * np.pi/width
    c_A[i] = 2*(tBcTop - tBcLeft) * (-(-1 + (-1)**(i+1))/(i+1)/np.pi)/np.sinh(lambda_A[i]*height)
    lambda_B[i] = (i+1) * np.pi/height
    c_B[i] = 2*(tBcRight - tBcBottom) * (-(-1 + (-1)**(i+1))/(i+1)/np.pi)/np.sinh(lambda_B[i]*width)

for i in range(np.size(x)):
    for j in range(np.size(y)):
        for k in range(numModes):
            T_A[j,i] += c_A[k] * np.sin(lambda_A[k] * x[i]) * np.sinh(lambda_A[k] * y[j])
            T_B[j,i] += c_B[k] * np.sin(lambda_B[k] * y[j]) * np.sinh(lambda_B[k] * x[i])

theta = T_A + T_B
T_analytic = theta + tBcBottom
TSliceHorizontalAnalytic = T_analytic[nY//2, :]
TSliceVerticalAnalytic = T_analytic[:, nX//2]

# relative errors for each center crossing solution line
deltaHorizontal = np.abs(TSliceHorizontal - TSliceHorizontalAnalytic)/TSliceHorizontalAnalytic
deltaVertical = np.abs(TSliceVertical - TSliceVerticalAnalytic)/TSliceVerticalAnalytic

# relative error plots
plt.figure(1)
plt.plot(x, deltaHorizontal)
plt.title("Relative Error in Horizontal Centerline")
plt.xlabel("X (m)")
plt.show()
plt.figure(2)
plt.plot(y, deltaVertical)
plt.title("Relative Error in Vertical Centerline")
plt.xlabel("Y (m)")
plt.show()


# plot the temperature distribution contour
plt.figure(3)
X = np.linspace(0, width, nX+2)
Y = np.linspace(0, height, nY+2)
X, Y = np.meshgrid(X, Y)
plt.contourf(X, Y, T_numerical, levels=50, cmap='plasma')
plt.colorbar(label = "Temperature (K)")
plt.title("Numerical Temperature Distribution")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()

# plot of analytic solution contour
plt.figure(4)
plt.contourf(X, Y, T_analytic, levels=50, cmap='plasma')
plt.colorbar(label = "Temperature (K)")
plt.title("Analytic Temperature Distribution")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.show()


