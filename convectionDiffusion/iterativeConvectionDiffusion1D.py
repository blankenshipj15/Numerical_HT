# Iterative solver for 1D convection diffusion problem with a leakage term using UDS and CDS schemes
# Solution provides inaccurate results for case (b) 

import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1.0
N = 20
dx = L / (N - 1)
x = np.linspace(0, L, N)
# dimensionless constants
Gamma = 1.0
A = 1.0
phi_left = 1.0
phi_right = 0.5

# Case a for inlet and outlet conditions
mx_out = 40 * Gamma * A / L
mx_in = 0.0
mL = (mx_out - mx_in) / L # assumed linear change in mx

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

# Function to calculate either the CDS or UDS for convection diffusion
def IterativeTdma(mx_in, mL, scheme="upwind", tol=1e-4, max_iter=1000):
    phi = np.linspace(phi_left, phi_right, N)  # initial guess for phi

    for iter in range(max_iter):
        a = np.zeros(N)
        b = np.zeros(N)
        c = np.zeros(N)
        d = np.zeros(N)

        # Interior nodes
        for i in range(1, N - 1):
            xP = x[i]
            xE = xP + dx / 2
            xW = xP - dx / 2

            m_e = mx_in - mL * xE
            m_w = mx_in - mL * xW

            D_e = Gamma * A / dx
            D_w = Gamma * A / dx

            F_e = m_e
            F_w = m_w

            if scheme == "upwind":
                c[i] = D_w + max(F_w, 0)
                b[i] = D_e + max(-F_e, 0)
                # Left boundary
                a[0] = Gamma*A/dx + mx_in - mL*(x[1]+dx/2)
                b[0] = Gamma*A/dx
                c[0] = 0.0
                d[0] = mx_in * phi_left
                # Right boundary
                a[-1] = Gamma*A/dx + mx_in - mL*(x[-1]+dx/2) 
                b[-1] = 0.0
                c[-1] = Gamma*A/dx + mx_in - mL*(x[-1]-dx/2)
                d[-1] = mx_out*phi_right
            elif scheme == "central":
                c[i] = D_w + 0.5 * F_w
                b[i] = D_e - 0.5 * F_e
                # left boundary
                a[0] = Gamma*A/dx + (mx_in - mL*(x[1]+dx/2))/2
                b[0] = Gamma*A/dx - (mx_in - mL*(x[1]+dx/2))
                c[0] = 0.0
                d[0] = mx_in * phi_left
                # right boundary
                a[-1] = Gamma*A/dx + mx_in - mL*(x[-1]+dx/2) 
                b[-1] = 0.0
                c[-1] = Gamma*A/dx + (mx_in - mL*(x[-1]-dx/2))
                d[-1] = mx_out*phi_right


            a[i] = b[i] + c[i] + (F_e - F_w)
            if mx_in == 0:
                d[i] = mL * phi_right * dx

            # the leakage term is dependent on phi at the current position in the domain 
            else:
                d[i] = mL * phi[i] * dx

        phi_new = tdmaSolver(a, b, c, d)

        if np.max(np.abs(phi_new - phi)) < tol:
            break
        phi = phi_new.copy()

    return phi

phiUpwind = IterativeTdma(mx_in, mL, scheme="upwind")
phiCentral = IterativeTdma(mx_in, mL, scheme="central")


plt.plot(x, phiUpwind, 'k-o', label='TDMA Upwind Case (b)')
plt.plot(x, phiCentral, 'r-o', label='TDMA Central Case (b)')
plt.xlabel('x')
plt.ylabel('$\\phi$')
plt.title('Case (b)')
plt.legend()
plt.show()
