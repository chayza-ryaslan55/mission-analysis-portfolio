
import numpy as np


def Eccentric_anomaly(M, e1, tol, max_iter):

    if not 0 <= e1 < 1:
        raise ValueError("Eccentricity must satisfy 0 <= e1 < 1.")

    if tol <= 0 or max_iter <= 0:
        raise ValueError("Tolerance and max_iter must be positive.")

    E = M

    for k in range(max_iter):
        f = E - e1 * np.sin(E) - M
        df = 1 - e1 * np.cos(E)

        if df == 0:
            raise RuntimeError("Derivative is zero; no convergence.")

        E_n = E - f / df

        # Residual evaluated at the NEW estimate
        residual = E_n - e1 * np.sin(E_n) - M

        print(f"Iteration {k + 1}: residual = {residual:.3e}")

        if abs(residual) <= tol:
            print(f"Converged in {k + 1} iterations")
            return E_n

        E = E_n

    raise RuntimeError("Newton-Raphson did not converge.")


def true_anomaly(E, e1):
    numerator = np.sqrt(1 + e1) * np.sin(E / 2)
    denominator = np.sqrt(1 - e1) * np.cos(E / 2)

    return 2 * np.arctan2(numerator, denominator)
