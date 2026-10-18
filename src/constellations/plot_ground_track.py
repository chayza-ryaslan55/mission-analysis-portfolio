import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add parent path to import propagation and geometry modules
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.abspath(os.path.join(current_dir, '..'))
if src_dir not in sys.path:
    sys.path.append(src_dir)
    
from propagation.propagation import propagate_orbit
from propagation.frames import keplerian_to_cartesian
from constellations.geometry import eci_to_ecef, ecef_to_geodetic, calculate_elevation


def plot_ground_track():
    print("Generating satellite trajectory for ground track visualization...")

    # 1. Define Orbital Elements (e.g., Sun-Synchronous Orbit)
    a = 7000.0  # Semi-major axis [km]
    e = 0.001  # Eccentricity
    i = np.radians(98.0)  # Inclination [rad] (SSO)
    raan = np.radians(0.0)  # RAAN [rad]
    arg_perigee = np.radians(0.0)  # Argument of perigee [rad]
    true_anomaly = np.radians(0.0)  # True anomaly [rad]

    initial_state = keplerian_to_cartesian(a, e, i, raan, arg_perigee, true_anomaly)

    # 2. Propagate for 3 full orbital periods to see multiple passes (~18,000 seconds)
    mu_earth = 398600.4418
    period = 2.0 * np.pi * np.sqrt((a ** 3) / mu_earth)
    t_max = 3.0 * period
    t_eval = np.linspace(0, t_max, 1500)

    sol = propagate_orbit(initial_state, t_span=(0, t_max), t_eval=t_eval)

    # 3. Extract Latitudes and Longitudes over time
    lats = []
    lons = []

    for idx, t in enumerate(sol.t):
        r_eci = sol.y[0:3, idx]
        r_ecef = eci_to_ecef(r_eci, t)
        lat, lon, alt = ecef_to_geodetic(r_ecef)
        lats.append(lat)
        lons.append(lon)

    lats = np.array(lats)
    lons = np.array(lons)

    # 4. Set up the 2D Matplotlib Figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot basic grid and map styling boundaries
    ax.set_xlim([-180, 180])
    ax.set_ylim([-90, 90])
    ax.set_xlabel("Longitude [degrees]", fontsize=12)
    ax.set_ylabel("Latitude [degrees]", fontsize=12)
    ax.set_title("Spacecraft Ground Track (3 Orbits)", fontsize=14, pad=15)
    ax.grid(True, linestyle='--', alpha=0.5)

    # Handle International Date Line wrapping to prevent horizontal cross-lines
    # We split the line segments whenever longitude jumps by more than 180 degrees
    lons_wrapped = lons.copy()
    lats_wrapped = lats.copy()

    # Plot trajectory segments
    segment_lons = [lons_wrapped[0]]
    segment_lats = [lats_wrapped[0]]

    for j in range(1, len(lons_wrapped)):
        if abs(lons_wrapped[j] - lons_wrapped[j - 1]) > 180.0:
            # Draw current segment and start a new one
            ax.plot(segment_lons, segment_lats, color='#1f77b4', linewidth=1.5)
            segment_lons = []
            segment_lats = []
        segment_lons.append(lons_wrapped[j])
        segment_lats.append(lats_wrapped[j])
    # Draw final segment
    ax.plot(segment_lons, segment_lats, color='#1f77b4', linewidth=1.5, label='Ground Track')

    # Mark start point and sub-satellite point highlights
    ax.scatter(lons[0], lats[0], color='green', s=60, zorder=5, label='Start Position')

    # Optional: Highlight a specific ground station (e.g., Toulouse)
    gs_lat, gs_lon = 43.6045, 1.4442
    ax.scatter(gs_lon, gs_lat, color='red', marker='^', s=80, zorder=5, label='Ground Station (Toulouse)')

    ax.legend(loc='lower left', frameon=True)
    plt.tight_layout()

    # Save output for README integration
    output_filename = 'ground_track_plot.png'
    plt.savefig(output_filename, dpi=300)
    print(f"Ground track plot saved successfully as '{output_filename}'!")

    plt.show()


if __name__ == "__main__":
    plot_ground_track()