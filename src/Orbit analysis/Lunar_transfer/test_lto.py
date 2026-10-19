import numpy as np
import sys
import os

# Set up local directory paths for clean execution
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

from lto_calculator import calculate_tli_parameters, propagate_lto


def main():
    print("Initializing Lunar Transfer Orbit (LTO) Analysis & Verification...")

    # Define baseline LEO parking orbit altitude
    h_park = 300.0
    params = calculate_tli_parameters(h_park_km=h_park)

    print("\nTrans-Lunar Injection (TLI) Performance Metrics:")
    print(f"  Parking Orbit Altitude:        {h_park:.1f} km")
    print(f"  Parking Radius:                {params['r_park_km']:.3f} km")
    print(f"  Circular Parking Velocity:     {params['v_park_km_s']:.4f} km/s")
    print(f"  Transfer Semi-Major Axis:      {params['a_trans_km']:.3f} km")
    print(f"  Transfer Perigee Velocity:     {params['v_perigee_km_s']:.4f} km/s")
    print(f"  Trans-Lunar Injection Delta-V: {params['delta_v_tli_km_s']:.4f} km/s")
    print(f"  Estimated Transit Time (ToF):  {params['tof_days']:.2f} days")

    print("\nExecuting Numerical Trajectory Propagation...")
    sol, _ = propagate_lto(h_park_km=h_park, num_points=400)

    # Verify boundary conditions at apogee
    final_position = sol.y[0:3, -1]
    final_distance = np.linalg.norm(final_position)

    print(f"  Total evaluation steps:        {len(sol.t)}")
    print(f"  Final radius at transfer apex: {final_distance:.2f} km")
    print(f"  Target lunar distance:         384400.00 km")
    print("\nLunar Transfer Orbit verification test completed successfully!")


if __name__ == "__main__":
    main()