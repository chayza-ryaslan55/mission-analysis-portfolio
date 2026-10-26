import numpy as np
from transfer import hohmann_transfer, lambert_universal_solver


def main():
    print("--- Testing Hohmann Transfer ---")
    r_initial = 6878.137  # 500 km altitude LEO [km]
    r_final = 42164.14  # Geostationary orbit (GEO) [km]

    hohmann = hohmann_transfer(r_initial, r_final)
    print(f"Transfer from {r_initial} km to {r_final} km:")
    print(f"  Delta-V 1:    {hohmann['dv1_km_s']:.4f} km/s")
    print(f"  Delta-V 2:    {hohmann['dv2_km_s']:.4f} km/s")
    print(f"  Total Delta-V: {hohmann['total_dv_km_s']:.4f} km/s")
    print(f"  Time of Flight: {hohmann['tof_hours']:.2f} hours")

    print("\n--- Testing Lambert Solver ---")
    # Define two positions and a time of flight
    r1 = np.array([7000.0, 0.0, 0.0])
    r2 = np.array([0.0, 7000.0, 1000.0])
    tof = 3600.0  # 1 hour transfer

    v1, v2 = lambert_universal_solver(r1, r2, tof)
    print(f"Initial Position: {r1}")
    print(f"Solved Initial Velocity (v1): {v1} km/s")
    print(f"Solved Final Velocity (v2):   {v2} km/s")
    print("\nLambert test completed successfully!")


if __name__ == "__main__":
    main()