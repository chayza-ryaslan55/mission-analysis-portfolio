import numpy as np
from scipy.integrate import solve_ivp

# Physical and orbital constants
MU_EARTH = 398600.4418  # Earth gravitational parameter [km^3/s^2]
R_EARTH = 6378.137  # Earth equatorial radius [km]
R_MOON = 384400.0  # Average Earth-Moon distance [km]


def calculate_tli_parameters(h_park_km=300.0, r_target_km=R_MOON):
    """
    Computes Trans-Lunar Injection (TLI) delta-V and transfer time of flight.

    Parameters:
    - h_park_km: Parking orbit altitude above Earth surface [km]
    - r_target_km: Target distance (Moon orbital radius) [km]

    Returns:
    - dict containing parking velocity, transfer semi-major axis, perigee velocity, TLI delta-V, and ToF
    """
    r_park = R_EARTH + h_park_km

    # 1. Circular parking orbit velocity
    v_park = np.sqrt(MU_EARTH / r_park)

    # 2. Transfer ellipse semi-major axis (perigee at parking radius, apogee at lunar radius)
    a_trans = (r_park + r_target_km) / 2.0

    # 3. Velocity at perigee of transfer ellipse using the Vis-viva equation
    v_perigee = np.sqrt(MU_EARTH * (2.0 / r_park - 1.0 / a_trans))

    # 4. Trans-Lunar Injection Delta-V requirement
    delta_v_tli = v_perigee - v_park

    # 5. Time of Flight (half of the elliptical transfer orbit period)
    tof_seconds = np.pi * np.sqrt((a_trans ** 3) / MU_EARTH)
    tof_days = tof_seconds / 86400.0

    return {
        "r_park_km": r_park,
        "v_park_km_s": v_park,
        "a_trans_km": a_trans,
        "v_perigee_km_s": v_perigee,
        "delta_v_tli_km_s": delta_v_tli,
        "tof_seconds": tof_seconds,
        "tof_days": tof_days
    }


def two_body_ode(t, state, mu):
    """
    Computes two-body equations of motion for numerical integration.

    Parameters:
    - t: Current time [s]
    - state: Combined position and velocity state vector [6x1]
    - mu: Gravitational parameter [km^3/s^2]

    Returns:
    - derivative state vector [6x1]
    """
    r = state[0:3]
    v = state[3:6]
    r_norm = np.linalg.norm(r)

    acceleration = -mu * r / (r_norm ** 3)
    return np.concatenate((v, acceleration))


def propagate_lto(h_park_km=300.0, num_points=1000):
    params = calculate_tli_parameters(h_park_km)
    r_park = params["r_park_km"]
    v_perigee = params["v_perigee_km_s"]

    # Initial state vector at perigee: start on +X, velocity pointing along -Y
    # so the ellipse sweeps across the positive X domain toward the Moon at +384,400 km
    r0 = np.array([r_park, 0.0, 0.0])
    v0 = np.array([0.0, -v_perigee, 0.0])
    initial_state = np.concatenate((r0, v0))

    t_span = (0.0, params["tof_seconds"])
    t_eval = np.linspace(0.0, params["tof_seconds"], num_points)

    sol = solve_ivp(
        two_body_ode,
        t_span,
        initial_state,
        args=(MU_EARTH,),
        t_eval=t_eval,
        rtol=1e-9,
        atol=1e-9
    )

    return sol, params