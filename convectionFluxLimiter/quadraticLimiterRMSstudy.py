# Problem: steady 1-D convection with source using FVM + flux limiter
import numpy as np
import matplotlib.pyplot as plt

# Parameters
L   = 10.0
u   = 1.0
rho = 1.0
A   = 10.0
phi_in = 100.0

# grid refinements (number of cells)
nX_list = [10, 20, 40, 80]

# ---------- TDMA solver ----------
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

# ---------- quadratic flux limiter ----------
def quadraticLimiter(r):
    # r = (phiE - phiP) / (phiP - phiW)
    if r <= 2:
        psi = (2*r + r**2) / (2 + r + r**2)
    else:
        psi = 1.0
    return psi

maxIter   = 1000
tolerance = 1e-5

dx_values   = []          # store dx for each grid
rms_errors  = []          # store RMS error for each grid
solutions   = {}          # store phi for each nX
x_centers   = {}          # store cell-center locations for each nX

for N in nX_list:
    # --- build grid for this refinement ---
    dx = L / N
    x  = np.linspace(dx/2, L - dx/2, N)   # cell centers
    phi_prev = np.ones(N)                 # initial guess

    # --- nonlinear iteration because of deferred correction ---
    for it in range(maxIter):
        a = np.zeros(N)
        b = np.zeros(N)
        c = np.zeros(N)
        d = np.zeros(N)

        for i in range(N):

            # exact cell-integrated source for the PDE:
            # d/dx(rho u phi) = (-2πA/L) sin(2πx/L)
            S_cell = A * (
                np.cos(2*np.pi*(x[i] + dx/2)/L) -
                np.cos(2*np.pi*(x[i] - dx/2)/L)
            )

            if i == 0:
                # inlet cell: upwind with phi(0) = phi_in
                a[i] = rho * u
                b[i] = 0.0
                c[i] = 0.0
                d[i] = rho * u * phi_in + S_cell

            elif i == 1:
                # second cell: still first-order upwind (no limiter yet)
                a[i] = rho * u
                b[i] = 0.0
                c[i] = rho * u
                d[i] = S_cell

            elif i < N - 1:
                # interior cells where we can use the limiter (need i-2, i-1, i, i+1)
                phiE_prev  = phi_prev[i+1]
                phiP_prev  = phi_prev[i]
                phiW_prev  = phi_prev[i-1]
                phiWW_prev = phi_prev[i-2]

                # gradient ratios with small epsilon to avoid division by zero
                rw = (phiP_prev - phiW_prev) / (phiW_prev - phiWW_prev + 1e-10)
                re = (phiE_prev - phiP_prev) / (phiP_prev - phiW_prev + 1e-10)

                psie = quadraticLimiter(re)
                psiw = quadraticLimiter(rw)

                # deferred correction = (high-order - upwind) flux contribution
                deferred_correction = (
                    rho*u * (phiW_prev + psiw * (phiW_prev - phiWW_prev) / 2.0 - phiW_prev)
                    - rho*u * (phiP_prev + psie * (phiP_prev - phiW_prev) / 2.0  - phiP_prev)
                )

                a[i] = rho * u
                b[i] = 0.0
                c[i] = rho * u
                d[i] = S_cell + deferred_correction

            else:
                # last cell: simple upwind outflow
                a[i] = rho * u
                b[i] = 0.0
                c[i] = rho * u
                d[i] = S_cell

        # solve for phi on this grid
        phi = tdmaSolver(a, b, c, d)

        # convergence check for deferred correction iterations
        residual = np.linalg.norm(phi - phi_prev, 2)
        if residual < tolerance:
            # print(f"N={N}, converged in {it} iterations, residual={residual:.3e}")
            break

        phi_prev = phi.copy()

    # store solution and grid
    solutions[N] = phi.copy()
    x_centers[N] = x.copy()

    # exact solution at cell centers:
    # phi(x) = 100 + (A/(rho u)) [cos(2πx/L) - 1]
    phi_exact_centers = phi_in + (A/(rho*u)) * (np.cos(2*np.pi*x/L) - 1.0)

    # RMS error at this grid
    rms = np.sqrt(np.mean((phi - phi_exact_centers)**2))
    dx_values.append(dx)
    rms_errors.append(rms)

# ---------- Plot solution for finest grid ----------
N_fine = nX_list[-1]
x_fine = x_centers[N_fine]
phi_fine = solutions[N_fine]

x_exact = np.linspace(0.0, L, 200)
phi_exact = phi_in + (A/(rho*u)) * (np.cos(2*np.pi*x_exact/L) - 1.0)

plt.figure()
plt.plot(np.concatenate(([0.0], x_fine)),
         np.concatenate(([phi_in], phi_fine)),
         'o-', label=f'Numerical (N={N_fine})')
plt.plot(x_exact, phi_exact, 'k--', label='Exact')
plt.xlabel('x')
plt.ylabel('phi')
plt.legend()
plt.grid(True)

# ---------- Plot RMS error vs dx ----------
dx_values = np.array(dx_values)
rms_errors = np.array(rms_errors)

plt.figure()
plt.loglog(dx_values, rms_errors, 'o-', label='RMS error')

# optional: reference first-order line through finest-grid point
ref_line = rms_errors[-1] * (dx_values / dx_values[-1])**1
plt.loglog(dx_values, ref_line, '--', label='O(dx) reference')

plt.gca().invert_xaxis()  # finer grids to the right if you like
plt.xlabel('dx')
plt.ylabel('RMS error')
plt.legend()
plt.grid(True, which='both')
plt.show()
