# coupledNonDimensional2D is a solution implementaion for a 2D transient heat conduction problem with Dirichlet BCs and TDMA solver.
# The problem considers coupled transient heat conduction between a metallic foam and a paraffin.
# The problem is discretized in time using a fully implicit scheme.
# Author: Jesse Blankenship
# Last Updated: 10/01/2025

import numpy as np
import matplotlib.pyplot as plt
from TDMA2DUnsteady import TDMA2DUnsteady
import time

# Problem parameters
H = 0.0 # interface parameter
epsilon = 0.8 # porosity
kRatio = 100.0 # kM/kP, thermal conductivity ratio
gammaM = 1.0 # non-dimensional thermal constant of metal foam
gammaP = 1.0/100.0    # alphaP / alphaM,  non-dimnesional thermal constant for paraffin
sCM = H/(1-epsilon) # constant source term for metal
sPM = H/(1-epsilon) # constant source term for metal
sCP = (H*kRatio*gammaP)/(epsilon) # constant source term for paraffin
sPP = (H*kRatio*gammaP)/(epsilon) # constant source term for paraffin
thetaBcLeft = 0.0 # Dirichlet left BC
thetaBcRight = 1.0  # Dirichlet right BC
thetaBcTop = 0.0 # Dirichlet top BC
thetaBcBottom = 1.0 # Dirichlet bottom BC
height = 1.0  # domain height (corresponds to non-dimensional y)
width = 1.0  # domain width (corresponds to non-dimensional x)
nX = 50 # number of cell centroids in x
nY = 50 # number of cell centroids in y
deltaX = width / (nX - 1) # distance between centroids in x
deltaY = height / (nY - 1) # distance between centroids in y
maxIter = 1000
tolerance = 1e-3 # convergence tolerance of solution at each timestep
timeStep = 0.01 # non-dimensional timestep size
steadyStateTolerance = 1e-6 # tolerance for steady state convergence in both materials
maxTime = 100
numberSteps = int(maxTime /  timeStep)

# Initialize non-dimensional temperature array and coefficient arrays
thetaM = np.zeros((nY, nX, numberSteps+1))
thetaP = np.zeros((nY, nX, numberSteps+1))
thetaPrevM = np.zeros((nY, nX)) # initial condition is 0 throughout domain
thetaPrevP = np.zeros((nY, nX)) # initial temperature distribution for paraffin
thetaNumericalM = np.zeros((nY+2, nX+2, numberSteps+1))
thetaNumericalP = np.zeros((nY+2, nX+2, numberSteps+1))

time_start = time.time()
t = 0.0 # initial time
step = 0
# need to begin with outer time loop
while t < maxTime:
    t += timeStep
    step += 1
    print(f"Step : {step}, Time : {t:.4f}")
    # first call the metal solver, initialized with initial conditions and assumed paraffin temp distribution
    tM = TDMA2DUnsteady(nX, nY, width, height, thetaPrevM, timeStep, gammaM, tolerance, maxIter, thetaBcLeft, thetaBcRight, thetaBcTop, thetaBcBottom, sCM, thetaPrevP, -sPM)
    # Solve using the solve method
    thetaM[:,:,step] = tM.solve()

    # at this point in time step, we have complete metal temperature distribution
    # now call paraffin solver using data from metal solver
    tP = TDMA2DUnsteady(nX, nY, width, height, thetaPrevP, timeStep, gammaP, tolerance, maxIter, thetaBcLeft, thetaBcRight, thetaBcTop, thetaBcBottom, sCP, thetaM[:,:,step], -sPP)
    thetaP[:,:,step] = tP.solve()

    # now we have a solution for both metal and paraffin at this time step
    # check for steady state convergences
    if np.max(np.abs(thetaP[:,:,step] - thetaPrevP)) < steadyStateTolerance and np.max(np.abs(thetaM[:,:,step] - thetaPrevM)) < steadyStateTolerance:
        print(f"Steady state reached at time {t}")
        break

    # append the BCs the solution for plotting
    thetaNumericalM[0,:,step] = thetaBcBottom
    thetaNumericalM[-1,:,step] = thetaBcTop
    thetaNumericalM[:,0,step] = thetaBcLeft
    thetaNumericalM[:,-1,step] = thetaBcRight
    thetaNumericalM[1:-1,1:-1,step] = thetaM[:,:,step]

    thetaNumericalP[0,:,step] = thetaBcBottom
    thetaNumericalP[-1,:,step] = thetaBcTop
    thetaNumericalP[:,0,step] = thetaBcLeft
    thetaNumericalP[:,-1,step] = thetaBcRight
    thetaNumericalP[1:-1,1:-1,step] = thetaP[:,:,step]

    # update the previous time step solution
    thetaPrevM = thetaM[:,:,step]
    thetaPrevP = thetaP[:,:,step]

time_end = time.time()
print("Computation time to converge: ", time_end - time_start)

# plot the non-dimensional steady state temperature distribution contour for metal
plt.figure(1)
X = np.linspace(0, width, nX+2)
Y = np.linspace(0, height, nY+2)
X, Y = np.meshgrid(X, Y)
plt.contourf(X, Y, thetaNumericalM[:,:,step-1], levels=50, cmap='plasma')
plt.colorbar(label = "Non-Dimensional Temperature ($\\theta$)")
plt.title("Non-Dimensional Metal Temperature Distribution")
plt.xlabel("x/D")
plt.ylabel("y/D")
plt.show()

# plot the non-dimensional steady state temperature distribution contour for paraffin
plt.figure(2)
X = np.linspace(0, width, nX+2)
Y = np.linspace(0, height, nY+2)
X, Y = np.meshgrid(X, Y)
plt.contourf(X, Y, thetaNumericalP[:,:,step-1], levels=50, cmap='plasma')
plt.colorbar(label = "Non-Dimensional Temperature ($\\theta$)")
plt.title("Non-Dimensional Paraffin Temperature Distribution")
plt.xlabel("x/D")
plt.ylabel("y/D")
plt.show()

# plot for centerpoint temperature over time
thetaCenterM = thetaNumericalM[nY//2,nX//2,:]
plt.figure(3)
plt.plot(np.linspace(0,t,step), thetaCenterM[:step])
plt.title("Centerpoint Metal Non-Dimensional Temperature")
plt.xlabel("Non-Dimensional Time ($\\tau$)")
plt.ylabel("Non-Dimensional Temperature ($\\theta$)")
plt.show()

# plot for centerpoint temperature over time
thetaCenterP = thetaNumericalP[nY//2,nX//2,:]
plt.figure(4)
plt.plot(np.linspace(0,t,step), thetaCenterP[:step])
plt.title("Centerpoint Paraffin Non-Dimensional Temperature")
plt.xlabel("Non-Dimensional Time ($\\tau$)")
plt.ylabel("Non-Dimensional Temperature ($\\theta$)")
plt.show()


# plot transient temperature distribution at vertical centerline
# plt.figure(3)
# plt.plot(Y[:,0], thetaNumerical[:,nX//2, 10], marker='o', markersize=4)
# plt.plot(Y[:,0], thetaNumerical[:,nX//2, 100], marker='.', markersize=4)
# plt.plot(Y[:,0], thetaNumerical[:,nX//2, 300], marker='x', markersize=4)
# plt.plot(Y[:,0], thetaNumerical[:,nX//2, 500], marker='^', markersize=4)
# plt.plot(Y[:,0], thetaNumerical[:,nX//2, 700], marker='s', markersize=4)
# plt.plot(Y[:,0], thetaNumerical[:,nX//2, step-1], marker='d', markersize=4)
# plt.legend(["$\\tau$ = 0.001", "$\\tau$ = 0.01","$\\tau$ = 0.03", "$\\tau$ = 0.05", "$\\tau$ = 0.07", "Steady State"])
# plt.title("Transient Non-Dimensional Temperature Distribution")
# plt.xlabel("Y/D")
# plt.ylabel("Non-Dimensional Temperature ($\\theta$)")
# plt.show()