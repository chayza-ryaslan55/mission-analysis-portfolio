import numpy as np

MU_EARTH = 398600.4418  # Earth gravitational parameter [km^3/s^2]


def keplerian_to_cartesian(a, e, i, raan, arg_perigee, true_anomaly):
    """
    Converts Classical Orbital Elements to ECI Cartesian state vectors.

    Parameters:
    - a: Semi-major axis [km]
    - e: Eccentricity [-]
    - i: Inclination [rad]
    - raan: Right Ascension of the Ascending Node [rad]
    - arg_perigee: Argument of Perigee [rad]
    - true_anomaly: True Anomaly [rad]

    Returns:
    - state: numpy array [x, y, z, vx, vy, vz] in [km, km/s]
    """
    # 1. Position and velocity in perifocal (orbital) frame
    p = a * (1.0 - e ** 2)  # Semi-latus rectum
    r_pf = (p / (1.0 + e * np.cos(true_anomaly))) * np.array([
        np.cos(true_anomaly),
        np.sin(true_anomaly),
        0.0
    ])

    v_pf = np.sqrt(MU_EARTH / p) * np.array([
        -np.sin(true_anomaly),
        e + np.cos(true_anomaly),
        0.0
    ])

    # 2. Rotation matrices to transform from perifocal to ECI frame
    # Using 3-1-3 or standard 3-1-3 Euler angle rotation sequence (RAAN, inclination, arg_perigee)
    cos_O, sin_O = np.cos(raan), np.sin(raan)
    cos_i, sin_i = np.cos(i), np.sin(i)
    cos_w, sin_w = np.cos(arg_perigee), np.sin(arg_perigee)

    # Transformation matrix from perifocal to ECI
    # R = R_z(-RAAN) * R_x(-i) * R_z(-arg_perigee) transposed correctly
    Q_sx = np.array([
        [cos_O * cos_w - sin_O * sin_w * cos_i, -cos_O * sin_w - sin_O * cos_w * cos_i, sin_O * sin_i],
        [sin_O * cos_w + cos_O * sin_w * cos_i, -sin_O * sin_w + cos_O * cos_w * cos_i, -cos_O * sin_i],
        [sin_w * sin_i, cos_w * sin_i, cos_i]
    ])

    r_eci = Q_sx @ r_pf
    v_eci = Q_sx @ v_pf

    return np.concatenate((r_eci, v_eci))