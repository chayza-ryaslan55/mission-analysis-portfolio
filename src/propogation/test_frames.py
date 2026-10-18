import numpy as np
from frames import keplerian_to_cartesian
from propagation import propagate_orbit


def main():
    print("Initializing Keplerian-to-Cartesian & Propagation Test...")

    # Define Classical Orbital Elements (COEs) for a Sun-Synchronous LEO
    a = 6878.137  # Semi-major axis [km] (approx 500 km altitude)
    e = 0.001  # Eccentricity (nearly circular) [-]
    i = np.radians(97.4)  # Inclination [rad] (Typical SSO inclination)
    raan = np.radians(30.0)  # Right Ascension of the Ascending Node [rad]
    arg_perigee = np.radians(0.0)  # Argument of Perigee [rad]
    true_anomaly = np.radians(0.0)  # True Anomaly [rad]

    # 1. Convert Classical Orbital Elements to Cartesian ECI State Vector
    initial_state = keplerian_to_cartesian(a, e, i, raan, arg_perigee, true_anomaly)

    print(f"\nConverted State Vector [km, km/s]:")
    print(f"  Position (r): {initial_state[0:3]}")
    print(f"  Velocity (v): {initial_state[3:6]}")

    # 2. Propagate the orbit for 1 orbital period (~5400 seconds)
    t_max = 5400.0
    t_eval = np.linspace(0, t_max, 100)

    print("\nRunning J2 Propagator...")
    sol = propagate_orbit(initial_state, t_span=(0, t_max), t_eval=t_eval)

    # 3. Display Results
    print(f"Propagation Status: {sol.message}")
    print(f"Initial Position [km]: {initial_state[0:3]}")
    print(f"Final Position [km]:   {sol.y[0:3, -1]}")
    print("\nTest completed successfully!")


if __name__ == "__main__":
    main()