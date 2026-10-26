import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from lto_calculator import propagate_lto, R_EARTH, R_MOON


def plot_lunar_transfer_3d(h_park=300.0, earth_scale=4.0, full_ellipse=True):
    sol, params = propagate_lto(h_park_km=h_park, num_points=5000)
    x, y, z = sol.y[:3]
    print(f"End position [km]: {sol.y[:3, -1]}")  # expect x ~ +384,400, y ~ 0

    r_p, r_a = params["r_park_km"], R_MOON
    a = params["a_trans_km"]
    e = (r_a - r_p) / (r_a + r_p)
    b = a * np.sqrt(1 - e**2)

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection="3d")

    # Earth
    r_draw = R_EARTH * earth_scale
    u, v = np.mgrid[0:2 * np.pi:40j, 0:np.pi:20j]
    ax.plot_surface(r_draw * np.cos(u) * np.sin(v),
                    r_draw * np.sin(u) * np.sin(v),
                    r_draw * np.cos(v), color="#1f77b4", alpha=0.8, linewidth=0)

    # Moon's orbit (reference circle)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(R_MOON * np.cos(th), R_MOON * np.sin(th), 0 * th,
            "--", color="gray", lw=0.8)

    # Full ellipse (dashed) so the whole orbit shape is visible
    if full_ellipse:
        t = np.linspace(0, 2 * np.pi, 600)
        c = a - r_p  # ellipse centre on the X axis
        ax.plot(c + a * np.cos(t), b * np.sin(t), 0 * t,
                ":", color="#ff7f0e", lw=1.2)

    # Simulated transfer arc and Moon at arrival
    ax.plot(x, y, z, color="#ff7f0e", lw=2.5)
    ax.scatter(x[-1], y[-1], z[-1], color="#7f7f7f", s=200, depthshade=False)

    # Limits from the data, with equal scaling on all axes
    pad = 30_000
    xlim = (-r_p - pad, r_a + pad)
    ylim = (-R_MOON * 0.35, R_MOON * 0.35)
    zlim = ylim
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_zlim(*zlim)
    ax.set_box_aspect([np.ptp(xlim), np.ptp(ylim), np.ptp(zlim)])
    ax.view_init(elev=35, azim=-60)

    ax.set_xlabel("X [km]", labelpad=12)
    ax.set_ylabel("Y [km]", labelpad=12)
    ax.set_zlabel("Z [km]", labelpad=12)
    ax.set_title(f"Lunar Transfer Orbit (ToF = {params['tof_days']:.2f} d, "
                 f"ΔV = {params['delta_v_tli_km_s']:.3f} km/s)")

    ax.legend(handles=[
        Line2D([0], [0], color="#ff7f0e", lw=2.5, label="Transfer arc (simulated)"),
        Line2D([0], [0], color="#ff7f0e", ls=":", label="Rest of ellipse"),
        Line2D([0], [0], color="gray", ls="--", label="Moon orbit"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#1f77b4",
               markersize=10, label=f"Earth (radius x{earth_scale:g})"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#7f7f7f",
               markersize=10, label="Moon (arrival point)"),
    ], loc="upper left")

    plt.tight_layout()
    plt.savefig("lto_trajectory_3d_plot.png", dpi=300)
    plt.show()


if __name__ == "__main__":
    plot_lunar_transfer_3d(earth_scale=2.0)

