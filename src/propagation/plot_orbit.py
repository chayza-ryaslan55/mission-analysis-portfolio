import numpy as np
import matplotlib.pyplot as plt
from frames import keplerian_to_cartesian
from propagation import R_EARTH, propagate_orbit


def plot_3d_orbit():
    print("Generating orbit trajectory for visualization...")

    # 1. Define Orbital Elements (e.g., an inclined LEO orbit)
    a = 7000.0  # Semi-major axis [km] (~622 km altitude)
    e = 0.001  # Eccentricity
    i = np.radians(45.0)  # Inclination [rad]
    raan = np.radians(10.0)  # RAAN [rad]
    arg_perigee = np.radians(0.0)  # Argument of perigee [rad]
    true_anomaly = np.radians(0.0)  # True anomaly [rad]

    # Convert to initial state vector
    initial_state = keplerian_to_cartesian(a, e, i, raan, arg_perigee, true_anomaly)

    # 2. Propagate for 2 full orbital periods (~10,800 seconds)
    # Estimate orbital period T = 2 * pi * sqrt(a^3 / mu)
    mu_earth = 398600.4418
    period = 2.0 * np.pi * np.sqrt((a ** 3) / mu_earth)
    t_max = 2.0 * period
    t_eval = np.linspace(0, t_max, 1000)

    sol = propagate_orbit(initial_state, t_span=(0, t_max), t_eval=t_eval)
    x, y, z = sol.y[0], sol.y[1], sol.y[2]

    # 3. Set up the 3D Matplotlib figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(projection='3d')

    # Plot the spacecraft trajectory
    ax.plot(x, y, z, label='Spacecraft Trajectory ($J_2$ Perturbed)', color='#1f77b4', linewidth=2)
    ax.scatter(x[0], y[0], z[0], color='green', s=50, label='Start Point')

    # Create a wireframe sphere representing the Earth
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_earth = R_EARTH * np.outer(np.cos(u), np.sin(v))
    y_earth = R_EARTH * np.outer(np.sin(u), np.sin(v))
    z_earth = R_EARTH * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_wireframe(x_earth, y_earth, z_earth, color='royalblue', alpha=0.15, rstride=4, cstride=4, label='Earth')

    # Formatting the plot
    ax.set_title("Spacecraft Orbit Propagation with $J_2$ Perturbation", fontsize=14, pad=20)
    ax.set_xlabel("X (ECI) [km]", labelpad=10)
    ax.set_ylabel("Y (ECI) [km]", labelpad=10)
    ax.set_zlabel("Z (ECI) [km]", labelpad=10)

    # Set equal aspect ratio for realistic orbital geometry
    max_range = max(np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z)))
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])

    ax.legend(loc='upper left')
    plt.tight_layout()

    # Save the figure so it can be added to your README later
    plt.savefig('orbit_3d_plot.png', dpi=300)
    print("Plot saved successfully as 'orbit_3d_plot.png'!")

    plt.show()


if __name__ == "__main__":
    plot_3d_orbit()