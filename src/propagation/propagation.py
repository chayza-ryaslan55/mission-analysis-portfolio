import numpy as np
from scipy.integrate import solve_ivp

# Physical constants for Earth (WGS-84 / standard astrodynamics constants)
MU_EARTH = 398600.4418  # Earth gravitational parameter [km^3/s^2]
R_EARTH = 6378.137      # Equatorial radius [km]
J2 = 1.08262668e-3      # J2 perturbation coefficient

def two_body_j2_dynamics(t, state):
    """
    Computes the equations of motion for a spacecraft in Earth orbit
    including the J2 perturbation.
    """
    r_vec = state[0:3]
    v_vec = state[3:6]

    r = np.linalg.norm(r_vec)
    x, y, z = r_vec

    # 1. Standard two-body acceleration
    a_two_body = -MU_EARTH * r_vec / (r ** 3)

    # 2. J2 Perturbation acceleration
    factor = (1.5 * J2 * MU_EARTH * (R_EARTH ** 2)) / (r ** 5)
    z_over_r_sq = (z / r) ** 2

    a_j2 = np.array([
        factor * x * (5.0 * z_over_r_sq - 1.0),
        factor * y * (5.0 * z_over_r_sq - 1.0),
        factor * z * (5.0 * z_over_r_sq - 3.0)
    ])

    # Total acceleration
    a_total = a_two_body + a_j2

    # Return state derivative [vx, vy, vz, ax, ay, az]
    return np.concatenate((v_vec, a_total))

def propagate_orbit(initial_state, t_span, t_eval=None, method='RK45', rtol=1e-9, atol=1e-9):
    """
    Propagates an orbit numerically using SciPy's IVP solver.
    """
    solution = solve_ivp(
        fun=two_body_j2_dynamics,
        t_span=t_span,
        y0=initial_state,
        t_eval=t_eval,
        method=method,
        rtol=rtol,
        atol=atol
    )
    return solution