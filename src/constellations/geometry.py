import numpy as np

# Earth constants (WGS-84)
R_EARTH = 6378.137  # Equatorial radius [km]
FLATTENING = 1.0 / 298.257223563  # Earth flattening factor
OMEGA_EARTH = 7.292115e-5  # Earth rotation rate [rad/s]


def eci_to_ecef(r_eci, t_sec, gmst_0=0.0):
    """
    Transforms an ECI position vector to ECEF accounting for Earth rotation.

    Parameters:
    - r_eci: Position vector in ECI [km] (3,)
    - t_sec: Elapsed time from initial epoch [seconds]
    - gmst_0: Initial Greenwich Mean Sidereal Time at t=0 [radians]

    Returns:
    - r_ecef: Position vector in ECEF [km] (3,)
    """
    # Current Earth rotation angle
    theta = gmst_0 + OMEGA_EARTH * t_sec

    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    # Rotation matrix about Z axis
    R3 = np.array([
        [cos_t, sin_t, 0.0],
        [-sin_t, cos_t, 0.0],
        [0.0, 0.0, 1.0]
    ])

    return R3 @ r_eci


def ecef_to_geodetic(r_ecef):
    """
    Converts ECEF Cartesian coordinates to geodetic latitude, longitude, and altitude (WGS-84).

    Parameters:
    - r_ecef: Position vector in ECEF [km] (3,)

    Returns:
    - lat: Latitude [degrees]
    - lon: Longitude [degrees]
    - alt: Altitude above reference ellipsoid [km]
    """
    x, y, z = r_ecef
    r = np.sqrt(x ** 2 + y ** 2)

    # Longitude
    lon = np.degrees(np.arctan2(y, x))

    # First and second eccentricities squared
    e2 = 2.0 * FLATTENING - FLATTENING ** 2

    # Iterative solution for latitude (Bowring's method)
    lat = np.arctan2(z, r * (1.0 - e2))  # initial guess
    for _ in range(5):
        sin_lat = np.sin(lat)
        N = R_EARTH / np.sqrt(1.0 - e2 * sin_lat ** 2)
        alt = r / np.cos(lat) - N
        lat = np.arctan2(z, r * (1.0 - e2 * (N / (N + alt))))

    lat = np.degrees(lat)
    return lat, lon, alt


def calculate_elevation(sat_ecef, gs_lat_deg, gs_lon_deg, gs_alt_km=0.0):
    """
    Calculates the elevation angle of a satellite relative to a ground station.

    Parameters:
    - sat_ecef: Satellite position in ECEF [km] (3,)
    - gs_lat_deg: Ground station latitude [degrees]
    - gs_lon_deg: Ground station longitude [degrees]
    - gs_alt_km: Ground station altitude [km]

    Returns:
    - elevation: Elevation angle [degrees]
    - azimuth: Azimuth angle [degrees]
    - range_km: Slant range distance [km]
    """
    lat_rad = np.radians(gs_lat_deg)
    lon_rad = np.radians(gs_lon_deg)

    # Convert ground station geodetic to ECEF
    e2 = 2.0 * FLATTENING - FLATTENING ** 2
    N = R_EARTH / np.sqrt(1.0 - e2 * np.sin(lat_rad) ** 2)

    gs_x = (N + gs_alt_km) * np.cos(lat_rad) * np.cos(lon_rad)
    gs_y = (N + gs_alt_km) * np.cos(lat_rad) * np.sin(lon_rad)
    gs_z = (N * (1.0 - e2) + gs_alt_km) * np.sin(lat_rad)
    gs_ecef = np.array([gs_x, gs_y, gs_z])

    # Relative position vector in ECEF
    rho_ecef = sat_ecef - gs_ecef

    # Transform ECEF relative vector to Topocentric horizon frame (ENU: East, North, Up)
    sin_lat, cos_lat = np.sin(lat_rad), np.cos(lat_rad)
    sin_lon, cos_lon = np.sin(lon_rad), np.cos(lon_rad)

    R_ecef_to_enu = np.array([
        [-sin_lon, cos_lon, 0.0],
        [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
        [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat]
    ])

    rho_enu = R_ecef_to_enu @ rho_ecef
    e, n, u = rho_enu

    # Slant range
    slant_range = np.linalg.norm(rho_enu)

    # Elevation and Azimuth
    elevation = np.degrees(np.arcsin(u / slant_range))
    azimuth = np.degrees(np.arctan2(e, n))
    if azimuth < 0:
        azimuth += 360.0

    return elevation, azimuth, slant_range