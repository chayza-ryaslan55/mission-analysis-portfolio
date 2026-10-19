import numpy as np

# Physical constants
BOLTZMANN_CONSTANT = 1.380649e-23  # Boltzmann constant [J/K]


def calculate_fspl(slant_range_km, frequency_ghz):
    """
    Computes Free Space Path Loss (FSPL) in decibels (dB).

    Parameters:
    - slant_range_km: Distance between satellite and ground station [km]
    - frequency_ghz: Carrier frequency [GHz]

    Returns:
    - fspl_db: Free space path loss [dB]
    """
    # Convert range to meters and frequency to Hz
    d_meters = slant_range_km * 1e3
    f_hz = frequency_ghz * 1e9

    # Speed of light
    c = 299792458.0

    # FSPL formula: (4 * pi * d * f / c)^2 expressed in dB
    fspl = (4.0 * np.pi * d_meters * f_hz / c) ** 2
    fspl_db = 10.0 * np.log10(fspl)

    return fspl_db


def compute_link_budget(
        power_tx_dbw,
        gain_tx_db,
        losses_tx_db,
        slant_range_km,
        frequency_ghz,
        gain_rx_db,
        losses_rx_db,
        system_temp_k,
        data_rate_kbps
):
    """
    Computes a complete downlink link budget.

    Parameters:
    - power_tx_dbw: Transmitter output power [dBW]
    - gain_tx_db: Transmitter antenna gain [dBi]
    - losses_tx_db: Transmitter implementation/feeder losses [dB]
    - slant_range_km: Slant range distance [km]
    - frequency_ghz: Carrier frequency [GHz]
    - gain_rx_db: Receiver ground station antenna gain [dBi]
    - losses_rx_db: Receiver atmospheric/feeder losses [dB]
    - system_temp_k: Receiver system noise temperature [K]
    - data_rate_kbps: Transmission data rate [kbps]

    Returns:
    - dict containing link parameters (EIRP, FSPL, Received Power, C/N0, Eb/N0)
    """
    # 1. Equivalent Isotropically Radiated Power (EIRP)
    eirp_dbw = power_tx_dbw + gain_tx_db - losses_tx_db

    # 2. Free Space Path Loss
    fspl_db = calculate_fspl(slant_range_km, frequency_ghz)

    # 3. Received Isotropic Power / Received Power (Pr)
    # Pr = EIRP - FSPL + G_rx - L_rx
    received_power_dbw = eirp_dbw - fspl_db + gain_rx_db - losses_rx_db

    # 4. Noise Power Spectral Density (N0) = k * T_sys [W/Hz]
    # Convert to dBW/Hz
    n0_dbw_hz = 10.0 * np.log10(BOLTZMANN_CONSTANT * system_temp_k)

    # 5. Carrier-to-Noise Density Ratio (C/N0) [dB-Hz]
    # C/N0 = Pr (in dBW) - N0 (in dBW/Hz)
    c_n0_dbhz = received_power_dbw - n0_dbw_hz

    # 6. Energy per Bit to Noise Spectral Density Ratio (Eb/N0) [dB]
    # Eb/N0 = C/N0 - 10 * log10(Data Rate in bps)
    data_rate_bps = data_rate_kbps * 1e3
    eb_n0_db = c_n0_dbhz - 10.0 * np.log10(data_rate_bps)

    return {
        "eirp_dbw": eirp_dbw,
        "fspl_db": fspl_db,
        "received_power_dbw": received_power_dbw,
        "c_n0_dbhz": c_n0_dbhz,
        "eb_n0_db": eb_n0_db
    }


def analyze_visibility_passes(sol, gs_lat, gs_lon, min_elevation=10.0):
    """
    Extracts continuous visibility windows (AOS and LOS) from propagated states.

    Parameters:
    - sol: SciPy OdeResult object from propagator
    - gs_lat, gs_lon: Ground station coordinates [degrees]
    - min_elevation: Minimum elevation angle mask [degrees]

    Returns:
    - passes: List of dicts detailing each pass start, end, max elevation, and duration.
    """
    from constellations.geometry import eci_to_ecef, ecef_to_geodetic, calculate_elevation

    passes = []
    in_pass = False
    current_pass = {}

    for idx, t in enumerate(sol.t):
        r_eci = sol.y[0:3, idx]
        r_ecef = eci_to_ecef(r_eci, t)
        el, az, rng = calculate_elevation(r_ecef, gs_lat, gs_lon)

        if el >= min_elevation:
            if not in_pass:
                # Acquisition of Signal (AOS)
                in_pass = True
                current_pass = {
                    "start_time_s": t,
                    "max_elevation": el,
                    "min_slant_range_km": rng
                }
            else:
                # Track maximum elevation during the pass
                if el > current_pass["max_elevation"]:
                    current_pass["max_elevation"] = el
                if rng < current_pass["min_slant_range_km"]:
                    current_pass["min_slant_range_km"] = rng
        else:
            if in_pass:
                # Loss of Signal (LOS)
                in_pass = False
                current_pass["end_time_s"] = t
                current_pass["duration_s"] = current_pass["end_time_s"] - current_pass["start_time_s"]
                passes.append(current_pass)

    # Handle case where simulation ends while still in a pass
    if in_pass:
        current_pass["end_time_s"] = sol.t[-1]
        current_pass["duration_s"] = current_pass["end_time_s"] - current_pass["start_time_s"]
        passes.append(current_pass)

    return passes