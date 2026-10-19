import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from lto_calculator import propagate_lto


def plot_lunar_transfer_3d():
    print("Generating 3D Lunar Transfer Orbit (LTO) visualization...")

    h_park = 300.0
    sol, params = propagate_lto(h_park_km=h_park, num_points=1000)

    r_earth = 6378.137
    r_moon = 384400.0

    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 1. Plot Earth at origin (using a wireframe sphere for 3D depth)
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 30)
    x_earth = r_earth * np.outer(np.cos(u), np.sin(v))
    y_earth = r_earth * np.outer(np.sin(u), np.sin(v))
    z_earth = r_earth * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x_earth, y_earth, z_earth, color='#1f77b4', alpha=0.6, label='Earth')

    # 2. Extract trajectory coordinates
    x = sol.y[0]
    y = sol.y[1]
    z = sol.y[2]

    # Plot Transfer Trajectory
    ax.plot(x, y, z, color='#ff7f0e', linewidth=2.0, label='Lunar Transfer Orbit (LTO)')

    # 3. Plot Moon position at apogee intercept
    ax.scatter([r_moon], [0.0], [0.0], color='#7f7f7f', s=150, zorder=5, label='Moon Position')

    # 4. Axis formatting and labels
    ax.set_xlabel("X Coordinate [km]", labelpad=10)
    ax.set_ylabel("Y Coordinate [km]", labelpad=10)
    ax.set_zlabel("Z Coordinate [km]", labelpad=10)
    ax.set_title("3D Trans-Lunar Injection (TLI) Trajectory Profile", fontsize=14, pad=20)

    # Set symmetric window bounds for proper spatial proportions
    ax.set_xlim([-50000, 420000])
    ax.set_ylim([-200000, 200000])
    ax.set_zlim([-200000, 200000])

    ax.legend(loc='upper left', frameon=True)
    plt.tight_layout()

    output_filename = 'lto_trajectory_3d_plot.png'
    plt.savefig(output_filename, dpi=300)
    print(f"3D Lunar transfer plot saved successfully as '{output_filename}'!")

    plt.show()


if __name__ == "__main__":
    plot_lunar_transfer_3d()