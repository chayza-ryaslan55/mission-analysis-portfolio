import numpy as np
import sys
import os

# Add the project root 'src' directory to the Python path dynamically
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '../..'))
if src_dir not in sys.path:
    sys.path.append(src_dir)

# Import modules using their package paths relative to 'src'
from propagation.propagation import propagate_orbit
from propagation.frames import keplerian_to_cartesian
from constellations.geometry import eci_to_ecef, ecef_to_geodetic, calculate_elevation


def main():
    print("Initializing Ground Track & Visibility Analysis Test...")

    # 1. Define LEO Orbit (e.g., Sun-synchronous orbit over Toulouse)
    a = 7000.0
    e = 0.001
    i = np.radians(98.0)
    raan = np.radians(0.0)
    arg_perigee = np.radians(0.0)
    true_anomaly = np.radians(0.0)

    initial_state = keplerian_to_cartesian(a, e, i, raan, arg_perigee, true_anomaly)

    # 2. Propagate for 1 full orbital period (~6000 seconds)
    t_max = 6000.0
    t_eval = np.linspace(0, t_max, 300)
    sol = propagate_orbit(initial_state, t_span=(0, t_max), t_eval=t_eval)

    # 3. Ground Station Definition (Toulouse, France)
    gs_lat = 43.6045
    gs_lon = 1.4442
    min_elevation = 10.0  # degrees mask angle

    visible_points = 0
    print(f"\nTracking satellite over ground station (Toulouse: {gs_lat}°N, {gs_lon}°E)...")

    for idx, t in enumerate(sol.t):
        r_eci = sol.y[0:3, idx]
        r_ecef = eci_to_ecef(r_eci, t)
        lat, lon, alt = ecef_to_geodetic(r_ecef)

        # Calculate visibility from ground station
        el, az, rng = calculate_elevation(r_ecef, gs_lat, gs_lon)

        if el >= min_elevation:
            visible_points += 1

    print(f"Total time steps evaluated: {len(sol.t)}")
    print(f"Steps where satellite is above {min_elevation}° elevation: {visible_points}")
    print("\nGeometry and Visibility test completed successfully!")


if __name__ == "__main__":
    main()