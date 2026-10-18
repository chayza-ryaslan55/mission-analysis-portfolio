import numpy as np

MU_EARTH = 398600.4418  # Earth gravitational parameter [km^3/s^2]


def hohmann_transfer(r1, r2):
    """
    Computes a coplanar Hohmann transfer between two circular orbits.

    Parameters:
    - r1: Radius of initial circular orbit [km]
    - r2: Radius of final circular orbit [km]

    Returns:
    - dict containing delta-v values and transfer time of flight.
    """
    # Circular orbital velocities
    v1 = np.sqrt(MU_EARTH / r1)
    v2 = np.sqrt(MU_EARTH / r2)

    # Semi-major axis of the transfer ellipse
    a_t = (r1 + r2) / 2.0

    # Velocities on the transfer ellipse at periapsis and apoapsis
    v_t1 = np.sqrt(MU_EARTH * (2.0 / r1 - 1.0 / a_t))
    v_t2 = np.sqrt(MU_EARTH * (2.0 / r2 - 1.0 / a_t))

    # Maneuver delta-Vs
    dv1 = abs(v_t1 - v1)
    dv2 = abs(v2 - v_t2)
    total_dv = dv1 + dv2

    # Time of flight (half of the elliptical period)
    tof = np.pi * np.sqrt((a_t ** 3) / MU_EARTH)

    return {
        "dv1_km_s": dv1,
        "dv2_km_s": dv2,
        "total_dv_km_s": total_dv,
        "tof_seconds": tof,
        "tof_hours": tof / 3600.0
    }


def lambert_universal_solver(r1_vec, r2_vec, tof, gm=MU_EARTH, prograde=True, max_iter=100, tol=1e-8):
    """
    Solves Lambert's boundary-value problem using the Universal Variable formulation.

    Parameters:
    - r1_vec: Initial position vector [km] (3,)
    - r2_vec: Final position vector [km] (3,)
    - tof: Time of flight between r1 and r2 [seconds]
    - gm: Gravitational parameter [km^3/s^2]
    - prograde: Boolean indicating direction (True for prograde, False for retrograde)

    Returns:
    - v1_vec: Initial velocity vector [km/s] (3,)
    - v2_vec: Final velocity vector [km/s] (3,)
    """
    r1 = np.linalg.norm(r1_vec)
    r2 = np.linalg.norm(r2_vec)

    cos_dtheta = np.dot(r1_vec, r2_vec) / (r1 * r2)
    cos_dtheta = np.clip(cos_dtheta, -1.0, 1.0)
    dtheta = np.arccos(cos_dtheta)

    # Determine direction of motion (short way vs long way)
    cross_z = r1_vec[0] * r2_vec[1] - r1_vec[1] * r2_vec[0]
    if prograde:
        if cross_z < 0:
            dtheta = 2.0 * np.pi - dtheta
    else:
        if cross_z >= 0:
            dtheta = 2.0 * np.pi - dtheta

    A = np.sin(dtheta) * np.sqrt((r1 * r2) / (1.0 - cos_dtheta))

    # Iterative solution for universal variable z
    # (Simplified robust Newton-Raphson loop for universal Kepler formulation)
    z = 0.0
    n_iter = 0

    while n_iter < max_iter:
        # Stumpff functions C(z) and S(z) approximations/evaluations
        if abs(z) > 1e-6:
            c = (1.0 - np.cos(np.sqrt(z))) / z
            s = (np.sqrt(z) - np.sin(np.sqrt(z))) / (z * np.sqrt(z))
        else:
            c = 0.5 - z / 24.0 + z ** 2 / 720.0
            s = (1.0 / 6.0) - z / 120.0 + z ** 2 / 5040.0

        y = r1 + r2 + A * (z * s - 1.0) / np.sqrt(c)

        if y < 0:
            # Shift z to avoid negative values during early iterations
            z += 1.0
            n_iter += 1
            continue

        x = np.sqrt(y / c)
        t_calc = (x ** 3 * s + A * np.sqrt(y)) / np.sqrt(gm)

        # Check convergence
        dt = t_calc - tof
        if abs(dt) < tol:
            break

        # Derivative dt/dz
        if abs(z) < 1e-6:
            dtdz = (np.sqrt(2.0) / 40.0) * (y ** 1.5) + (A / 8.0) * (np.sqrt(y) + A * np.sqrt(1.0 / (2.0 * y)))
        else:
            dtdz = (x ** 3 * (s - 3.5 * c * s / (2.0 * c) + 0.75 * z / c) + 0.5 * A * (y / c) ** 0.5) / (
                        z * gm)  # simplified gradient step

        z -= dt / max(abs(dtdz), 1e-3)
        n_iter += 1

    # Compute Lagrange coefficients f, g, f_dot, g_dot
    f = 1.0 - y / r1
    g = A * np.sqrt(y / gm)
    g_dot = 1.0 - y / r2

    v1_vec = (r2_vec - f * r1_vec) / g
    v2_vec = (g_dot * r2_vec - r1_vec) / g

    return v1_vec, v2_vec