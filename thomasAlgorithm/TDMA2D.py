# TDMA class used for solving 2D steady state heat conduction problems
# coordinates are assumed cartesian with uniform grid spacing
# boundary conditions are assumed Dirichlet
# this class uses line by line solving, sweeping both horizontally and vertically
# Author: Jesse Blankenship
# Last Updated: 9/21/2025

import numpy as np

class TDMA2D:

    def __init__ (self, nX, nY, width, height, gamma, tol, maxIter, TLeft, TRight, TTop, TBottom):
        self.nX = nX
        self.nY = nY
        self.deltaX = width/(nX-1)
        self.deltaY = height/(nY-1)
        self.gamma = gamma
        self.tol = tol
        self.maxIter = maxIter
        self.TLeft = TLeft
        self.TRight = TRight
        self.TTop = TTop
        self.TBottom = TBottom
        self.T_prev = np.zeros((nY, nX))  # initialize temperature array
        self.T = np.zeros((nY, nX))  # solution array

        # the coefficients are sized assuming a square domain
        self.a = np.zeros((nY, 1))
        self.b = np.zeros((nY, 1))
        self.c = np.zeros((nY, 1))
        self.d = np.zeros((nY, 1))

    def tdmaSolver(self, a, b, c, d):
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
    
    def horizontalSweep(self):
        nX = self.nX
        nY = self.nY
        deltaX = self.deltaX
        deltaY = self.deltaY
        gamma = self.gamma
        T_prev = self.T_prev
        T = self.T
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        tBcLeft = self.TLeft
        tBcRight = self.TRight
        tBcTop = self.TTop
        tBcBottom = self.TBottom
        for i in range(nX):
            if i == 0: # this is the first vertical line
                for j in range(nY): # build the coefficients for the first vertical line
                    if j == 0: # bottom left corner
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2) + gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2))
                        b[j] = gamma*deltaX/deltaY
                        c[j] = 0
                        d[j] = (gamma*deltaX/(deltaY/2))*tBcBottom + (gamma*deltaY/deltaX)*T_prev[j,i+1] +(gamma*deltaY/(deltaX/2))*tBcLeft
                    elif j == nY-1: # top left corner
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2) + gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2))
                        b[j] = 0
                        c[j] = gamma*deltaX/deltaY
                        d[j] = (gamma*deltaX/(deltaY/2))*tBcTop + (gamma*deltaY/deltaX) * T_prev[j,i+1] +(gamma*deltaY/(deltaX/2))*tBcLeft
                    else: # left edge
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/deltaY + gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2))
                        b[j] = gamma*deltaX/deltaY
                        c[j] = gamma*deltaX/deltaY
                        d[j] = (gamma*deltaY/deltaX)*T_prev[j,i+1] +(gamma*deltaY/(deltaX/2))*tBcLeft
                # solve the tridiagonal system for this vertical line
                T[:,i] = self.tdmaSolver(a.flatten(), b.flatten(), c.flatten(), d.flatten())
            elif i == nX-1: # this is the last vertical line
                for j in range(nY): # build the coefficients for the last vertical line
                    if j == 0: # bottom right corner
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2) + gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2))
                        b[j] = 0
                        c[j] = gamma*deltaX/deltaY
                        d[j] = (gamma*deltaX/(deltaY/2))*tBcBottom + (gamma*deltaY/deltaX)*T_prev[j,i-1] +(gamma*deltaY/(deltaX/2))*tBcRight
                    elif j == nY-1: # top right corner
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2) + gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2))
                        b[j] = gamma*deltaX/deltaY
                        c[j] = 0
                        d[j] = (gamma*deltaX/(deltaY/2))*tBcTop + (gamma*deltaY/deltaX)*T_prev[j,i-1] +(gamma*deltaY/(deltaX/2))*tBcRight
                    else: # right edge
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/deltaY + gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2))
                        b[j] = 0
                        c[j] = gamma*deltaX/deltaY
                        d[j] = (gamma*deltaY/deltaX)*T_prev[j,i-1] +(gamma*deltaY/(deltaX/2))*tBcRight
                # solve the tridiagonal system for this vertical line
                T[:,i] = self.tdmaSolver(a.flatten(), b.flatten(), c.flatten(), d.flatten())
            else: # this is an interior vertical line
                for j in range(nY): # build the coefficients for this interior vertical line
                    if j == 0: # bottom edge
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2) + gamma*deltaY/deltaX + gamma*deltaY/deltaX)
                        b[j] = gamma*deltaX/deltaY
                        c[j] = 0
                        d[j] = (gamma*deltaX/(deltaY/2))*tBcBottom + (gamma*deltaY/deltaX)*T_prev[j,i+1] + (gamma*deltaY/deltaX)*T_prev[j,i-1]
                    elif j == nY-1: # top edge
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2) + gamma*deltaY/deltaX + gamma*deltaY/deltaX)
                        b[j] = 0
                        c[j] = gamma*deltaX/deltaY
                        d[j] = (gamma*deltaX/(deltaY/2))*tBcTop + (gamma*deltaY/deltaX)*T_prev[j,i+1] + (gamma*deltaY/deltaX)*T_prev[j,i-1]
                    else: # interior points
                        a[j] = (gamma*deltaX/deltaY + gamma*deltaX/deltaY + gamma*deltaY/deltaX + gamma*deltaY/deltaX)
                        b[j] = gamma*deltaX/deltaY
                        c[j] = gamma*deltaX/deltaY
                        d[j] = (gamma*deltaY/deltaX)*T_prev[j,i+1] + (gamma*deltaY/deltaX)*T_prev[j,i-1]
                # solve the tridiagonal system for this vertical line
                T[:,i] = self.tdmaSolver(a.flatten(), b.flatten(), c.flatten(), d.flatten())
        return T

    def verticalSweep(self):
        nX = self.nX
        nY = self.nY
        deltaX = self.deltaX
        deltaY = self.deltaY
        gamma = self.gamma
        T = self.T
        a = self.a
        b = self.b
        c = self.c
        d = self.d
        tBcLeft = self.TLeft
        tBcRight = self.TRight
        tBcTop = self.TTop
        tBcBottom = self.TBottom
        for j in range(nY):
                if j == 0: # this is the first (bottom) horizontal line
                    for i in range(nX): # build the coefficients for the first horizontal line
                        if i == 0: # bottom left corner
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2) + gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2))
                            b[i] = gamma*deltaY/deltaX
                            c[i] = 0
                            d[i] = (gamma*deltaY/(deltaX/2))*tBcLeft + (gamma*deltaX/deltaY)*T[j+1,i] +(gamma*deltaX/(deltaY/2))*tBcBottom
                        elif i == nX-1: # bottom right corner
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2) + gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2))
                            b[i] = 0
                            c[i] = gamma*deltaY/deltaX
                            d[i] = (gamma*deltaY/(deltaX/2))*tBcRight + (gamma*deltaX/deltaY)*T[j+1,i] +(gamma*deltaX/(deltaY/2))*tBcBottom
                        else: # bottom edge
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/deltaX + gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2))
                            b[i] = gamma*deltaY/deltaX
                            c[i] = gamma*deltaY/deltaX
                            d[i] = (gamma*deltaX/deltaY)*T[j+1,i] + (gamma*deltaX/(deltaY/2))*tBcBottom
                    # solve the tridiagonal system for this horizontal line
                    T[j,:] = self.tdmaSolver(a.flatten(), b.flatten(), c.flatten(), d.flatten())
                elif j == nY-1: # this is the last horizontal line (top)
                    for i in range(nX): # build the coefficients for the last horizontal line
                        if i == 0: # top left corner
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2) + gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2))
                            b[i] = gamma*deltaY/deltaX
                            c[i] = 0
                            d[i] = (gamma*deltaY/(deltaX/2))*tBcLeft + (gamma*deltaX/deltaY)*T[j-1,i] +(gamma*deltaX/(deltaY/2))*tBcTop
                        elif i == nX-1: # top right corner
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2) + gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2))
                            b[i] = 0
                            c[i] = gamma*deltaY/deltaX
                            d[i] = (gamma*deltaY/(deltaX/2))*tBcRight + (gamma*deltaX/deltaY)*T[j-1,i] +(gamma*deltaX/(deltaY/2))*tBcTop
                        else: # top edge
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/deltaX + gamma*deltaX/deltaY + gamma*deltaX/(deltaY/2))
                            b[i] = gamma*deltaY/deltaX
                            c[i] = gamma*deltaY/deltaX
                            d[i] = (gamma*deltaX/deltaY)*T[j-1,i] +(gamma*deltaX/(deltaY/2))*tBcTop
                    # solve the tridiagonal system for this horizontal line
                    T[j,:] = self.tdmaSolver(a.flatten(), b.flatten(), c.flatten(), d.flatten())
                else: # this is an interior horizontal line
                    for i in range(nX): # build the coefficients for this interior horizontal line
                        if i == 0: # left edge
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2) + gamma*deltaX/deltaY + gamma*deltaX/deltaY)
                            b[i] = gamma*deltaY/deltaX
                            c[i] = 0
                            d[i] = (gamma*deltaY/(deltaX/2))*tBcLeft + (gamma*deltaX/deltaY)*T[j+1,i] + (gamma*deltaX/deltaY)*T[j-1,i]
                        elif i == nX-1: # right edge
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/(deltaX/2) + gamma*deltaX/deltaY + gamma*deltaX/deltaY)
                            b[i] = 0
                            c[i] = gamma*deltaY/deltaX
                            d[i] = (gamma*deltaY/(deltaX/2))*tBcRight + (gamma*deltaX/deltaY)*T[j+1,i] + (gamma*deltaX/deltaY)*T[j-1,i]
                        else: # interior points
                            a[i] = (gamma*deltaY/deltaX + gamma*deltaY/deltaX + gamma*deltaX/deltaY + gamma*deltaX/deltaY)
                            b[i] = gamma*deltaY/deltaX
                            c[i] = gamma*deltaY/deltaX
                            d[i] = (gamma*deltaX/deltaY)*T[j+1,i] + (gamma*deltaX/deltaY)*T[j-1,i]
                    # solve the tridiagonal system for this horizontal line
                    T[j,:] = self.tdmaSolver(a.flatten(), b.flatten(), c.flatten(), d.flatten())
        return T

    def solve(self):
        # method to sweep horizontally and vertically
        maxIter = self.maxIter
        iter = 0
        while iter < maxIter:
            self.horizontalSweep()
            self.verticalSweep()
            # convergence check
            if np.max(np.abs(self.T - self.T_prev)) < self.tol:
                break

            self.T_prev = np.copy(self.T)
            iter += 1

        return self.T
