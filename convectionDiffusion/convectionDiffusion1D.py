import numpy as np
import matplotlib.pyplot as plt

# ------------------------------
# Parameters
# ------------------------------
L = 1.0
N = 21
dx = L / (N - 1)
x = np.linspace(0, L, N)
Gamma = 1.0
A = 1.0
phi_left = 1.0
phi_right = 0.5
phi_leak = 0.0  # ambient/leak value

# ------------------------------
# TDMA Solver
# ------------------------------
def tdmaSolver(a, b, c, d):
    n = len(a)
    P = np.zeros(n)
    Q = np.zeros(n)
    phi = np.zeros(n)

    P[0] = b[0] / a[0]
    Q[0] = d[0] / a[0]

    for i in range(1, n):
        denom = a[i] - c[i] * P[i - 1]
        P[i] = b[i] / denom
        Q[i] = (d[i] + c[i] * Q[i - 1]) / denom

    phi[-1] = Q[-1]
    for i in range(n - 2, -1, -1):
        phi[i] = P[i] * phi[i + 1] + Q[i]

    return phi


# ------------------------------
# Iterative TDMA Solver
# ------------------------------
def IterativeTdma(mx_in, mx_out, scheme="upwind", tol=1e-6, max_iter=1000):
    mL = (mx_in - mx_out) / L  # slope of mx(x)
    phi = np.linspace(phi_left, phi_right, N)  # initial guess

    for it in range(max_iter):
        a = np.zeros(N)
        b = np.zeros(N)
        c = np.zeros(N)
        d = np.zeros(N)

        # Dirichlet BCs
        a[0], d[0] = 1.0, phi_left
        a[-1], d[-1] = 1.0, phi_right

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

            # Coefficients for UDS or CDS
            if scheme == "upwind":
                c[i] = D_w + max(F_w, 0)
                b[i] = D_e + max(-F_e, 0)
            elif scheme == "central":
                c[i] = D_w + 0.5 * F_w
                b[i] = D_e - 0.5 * F_e
            else:
                raise ValueError("scheme must be 'upwind' or 'central'")

            # Leakage term: acts as a sink
            a[i] = b[i] + c[i] + (F_e - F_w) + abs(mL) * A * dx
            d[i] = 0.0  # no fixed source term

        phi_new = tdmaSolver(a, b, c, d)

        if np.max(np.abs(phi_new - phi)) < tol:
            print(f"{scheme} converged in {it} iterations")
            break
        phi = phi_new.copy()

    return phi



# ------------------------------
# Case (a): mx_in = 40ΓA, mx_out = 0
# ------------------------------
mx_in_a = 40 * Gamma * A
mx_out_a = 0.0
phiUpwind_a = IterativeTdma(mx_in_a, mx_out_a, scheme="upwind")
phiCentral_a = IterativeTdma(mx_in_a, mx_out_a, scheme="central")

# ------------------------------
# Case (b): mx_in = 0, mx_out = 40ΓA
# ------------------------------
mx_in_b = 0.0
mx_out_b = 40 * Gamma * A
phiUpwind_b = IterativeTdma(mx_in_b, mx_out_b, scheme="upwind")
phiCentral_b = IterativeTdma(mx_in_b, mx_out_b, scheme="central")

# ------------------------------
# Plot results
# ------------------------------
plt.figure(figsize=(8,5))
plt.plot(x, phiUpwind_a, 'ko-', label='Upwind Case (a)')
plt.plot(x, phiCentral_a, 'r--o', label='Central Case (a)')
plt.plot(x, phiUpwind_b, 'bs--', label='Upwind Case (b)')
plt.plot(x, phiCentral_b, 'g-.s', label='Central Case (b)')
plt.xlabel('x')
plt.ylabel('$\\phi$')
plt.title('1D Convection–Diffusion with Leakage Term')
plt.legend()
plt.grid(True)
plt.show()
