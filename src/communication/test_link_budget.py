import numpy as np
import sys
import os

# Set up paths to import sibling packages
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '..'))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from propagation.propagation import propagate_orbit
from propagation.frames import keplerian_to_cartesian
from communication.link_budget import compute_link_budget, analyze_visibility_passes


def main():
    print("Initializing Link Budget and Pass Analysis Test...")

    # 1. Define LEO Orbit and propagate
    a = 7000.0
    e = 0.001
    i = np.radians(98.0)
    initial_state = keplerian_to_cartesian(a, e, i, 0.0, 0.0, 0.0)

    t_max = 6000.0
    t_eval = np.linspace(0, t_max, 600)
    sol = propagate_orbit(initial_state, t_span=(0, t_max), t_eval=t_eval)

    # 2. Ground Station (Toulouse)
    gs_lat = 43.6045
    gs_lon = 1.4442

    # Extract visibility passes
    passes = analyze_visibility_passes(sol, gs_lat, gs_lon, min_elevation=10.0)
    print(f"\nFound {len(passes)} visibility pass(es) over Toulouse during simulation.")

    if len(passes) > 0:
        p = passes[0]
        print(f"Pass Details:")
        print(f"  Duration: {p['duration_s']:.1f} seconds")
        print(f"  Max Elevation: {p['max_elevation']:.2f} degrees")
        print(f"  Minimum Slant Range: {p['min_slant_range_km']:.2f} km")

        # 3. Compute Link Budget at Closest Approach / Minimum Slant Range
        freq_ghz = 2.225  # S-band downlink frequency [GHz]
        budget = compute_link_budget(
            power_tx_dbw=3.0,  # 2 Watts transmitter power (~5 dBW)
            gain_tx_db=2.15,  # Dipole antenna gain [dBi]
            losses_tx_db=1.0,  # Feeder/cable losses [dB]
            slant_range_km=p['min_slant_range_km'],
            frequency_ghz=freq_ghz,
            gain_rx_db=35.0,  # Ground station dish gain [dBi]
            losses_rx_db=0.5,  # Atmospheric/polarization losses [dB]
            system_temp_k=150.0,  # System noise temperature [K]
            data_rate_kbps=100.0  # Downlink rate [kbps]
        )

        print("\nDownlink Link Budget Results (at max elevation):")
        print(f"  Free Space Path Loss (FSPL): {budget['fspl_db']:.2f} dB")
        print(f"  Received Power:              {budget['received_power_dbw']:.2f} dBW")
        print(f"  Carrier-to-Noise (C/N0):     {budget['c_n0_dbhz']:.2f} dB-Hz")
        print(f"  Energy per Bit (Eb/N0):      {budget['eb_n0_db']:.2f} dB")

    print("\nLink budget and visibility test completed successfully!")


if __name__ == "__main__":
    main()