# SIMPLE algorithm to calculate velocity and pressure field for a 1D steady nozzle
# the momenetum equations is linearized
import numpy as np

# parameters
rho = 1.0
AA = 3.0
AB = 1.0
P1 = 28.0
P3 = 0.0
alphaU = 0.8
alphaP = 0.8

# Inital guesses
uAprev = 5.0/3.0
uBprev = 5.0
P2prev = 25.0

maxIter = 1000
for iter in range(maxIter):
    # solve momentum to get starred velocities (step 2 of SIMPLE)
    dA = (alphaU*AA)/(rho*AA*uAprev)
    bA = uAprev*(1-alphaU)/(rho*AA*uAprev)
    uAstar = dA*(P1 - P2prev) + bA

    dB = (alphaU*AB)/(rho*AB*uBprev)
    bB = uBprev*(1-alphaU)/(rho*AB*uBprev)
    uBstar = dB*(P2prev - P3) + bB

    # solve pressure correction
    P2prime = (uAstar*AA - uBstar*AB)/(dB*AB + dA*AA)
    P2new = P2prev + alphaP*P2prime

    # update velocity
    uAnew = uAstar - dA*P2prime
    uBnew = uBstar + dB*P2prime

    # check continutity
    # continuity = abs(uAnew*AA - uBnew*AB)
    # if continuity < 1e-6:
    #     print(f'Converged in {iter + 1} iterations. uA = {uAnew:.4f}, uB = {uBnew:.4f}, P2 = {P2new:.4f}')
    #     break

    # update for next iteration
    uAprev = uAnew
    uBprev = uBnew
    P2prev = P2new

print(f'uA = {uAnew:.4f}, uB = {uBnew:.4f}, P2 = {P2new:.4f}')

#####################################################
## Solution is uA = 1.5304, uB = 4.5911, P2 = 25.2 ##
#####################################################