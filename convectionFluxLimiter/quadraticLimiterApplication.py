# Problem solving steady 1-D convection equation with a source using finite volume method
# Using quadratic flux limiter for higher order accuracy as compared to first order upwind
import numpy as np
import matplotlib.pyplot as plt

# Parameters
L   = 10.0
u   = 1.0
rho = 1.0
A   = 10.0         
nX  = 50
dx  = L / (nX - 1) 
x   = np.linspace(dx/2, L - dx/2, nX) # locating the cell centers
nRefinements = 5

phi_prev = np.ones(nX)  # initial guess

# TDMA solver
def tdmaSolver(a, b, c, d):
    nPoints = np.size(a)

    P = np.zeros(nPoints)
    Q = np.zeros(nPoints)
    phi = np.zeros(nPoints)

    # forward sweep
    P[0] = b[0] / a[0]
    Q[0] = d[0] / a[0]
    for i in range(1, nPoints):
        denom = a[i] - c[i] * P[i-1]
        P[i]  = b[i] / denom
        Q[i]  = (d[i] + c[i] * Q[i-1]) / denom

    # back substitution
    phi[-1] = Q[-1]
    for i in range(nPoints-2, -1, -1):
        phi[i] = P[i] * phi[i+1] + Q[i]

    return phi

# Quadratic flux limiter function
def quadraticLimiter(r):
    # r = (phiE - phiP) / (phiP - phiW)
    if r <= 2:
        psi = (2*r + r**2) / (2 + r + r**2)
    else:
        psi = 1.0
    return psi

maxIter   = 1000
# convergence tolerance for TDMA iterations
tolerance = 1e-5
for iter in range(maxIter):
    # allocate tri-diagonal coefficients each iteration
    a = np.zeros(nX)
    b = np.zeros(nX)   
    c = np.zeros(nX)   
    d = np.zeros(nX)   

    for i in range(nX):
        
        S_cell = A * (np.cos(2*np.pi*(x[i] + dx/2)/L) - np.cos(2*np.pi*(x[i] - dx/2)/L))

        if i == 0:
            # first order upwind with inlet bc
            a[i] = rho * u
            b[i] = 0.0
            c[i] = 0.0
            d[i] = rho * u * 100.0 + S_cell

        elif i == 1:
            # first order upwind still
            a[i] = rho * u      
            b[i] = 0.0
            c[i] = rho * u      
            d[i] = S_cell

        elif i < nX - 1:
            # interior cells where we can use the limiter (i-2, i-1, i, i+1)
            phiE_prev  = phi_prev[i+1]
            phiP_prev  = phi_prev[i]
            phiW_prev  = phi_prev[i-1]
            phiWW_prev = phi_prev[i-2]

            rw = (phiP_prev - phiW_prev) / (phiW_prev - phiWW_prev + 1e-10)
            re = (phiE_prev - phiP_prev) / (phiP_prev - phiW_prev + 1e-10)

            psie = quadraticLimiter(re)
            psiw = quadraticLimiter(rw)

            # deferred correction using quadratic limiter
            deferred_correction = (rho*u * (phiW_prev + psiw * (phiW_prev - phiWW_prev) / 2.0 - phiW_prev) 
            - rho*u * (phiP_prev + psie * (phiP_prev - phiW_prev) / 2.0  - phiP_prev))

            a[i] = rho * u
            b[i] = 0.0
            c[i] = rho * u
            d[i] = S_cell + deferred_correction

        else:
            a[i] = rho * u
            b[i] = 0.0
            c[i] = rho * u
            d[i] = S_cell

    phi = tdmaSolver(a, b, c, d)

    # convergence check
    residual = np.linalg.norm(phi - phi_prev, 2)
    
    if residual < tolerance:
        break

    phi_prev = phi.copy()

# Plotting results
x_plotting = np.linspace(0, L, nX+1)

phi_exact = np.zeros_like(x_plotting)
phi_exact[:] = 100.0 + (A/(rho*u)) * (np.cos(2*np.pi*x_plotting/L)-1)

phi_plotting = np.zeros(nX+1)
phi_plotting[0] = 100.0
phi_plotting[1:] = phi

plt.plot(x_plotting, phi_plotting, marker='o',
         label='Numerical (N = 50)')
plt.plot(x_plotting, phi_exact, 'r--', label='Exact Solution')
plt.xlabel('x')
plt.ylabel('phi')
plt.legend()
plt.grid(True)
plt.show()

# also plotting RMS error between numerical and exact solution to identify effect of grid refinement
rms_error = np.sqrt(np.mean((phi_plotting - phi_exact)**2))
print(f"RMS Error: {rms_error:.6f}")
rms_errors = np.zeros(nRefinements-1)
refinements = np.array([10, 50, 100, 200])
rms_errors[0] = 3.970431
rms_errors[1] = 0.480292
rms_errors[2] = 0.222224
rms_errors[3] = 0.107079

plt.plot(refinements,rms_errors)
plt.xlabel('# Grid Points (N)')
plt.ylabel('RMS Error')
plt.show()