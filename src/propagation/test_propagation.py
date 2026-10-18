import numpy as np

# Import from the local propagator module
from propagation import propagate_orbit

def main():
    print("Initializing J2 Orbit Propagation Test...")

    # Example LEO state vector [km, km/s] (approx. 500 km altitude, inclined orbit)
    r_init = [6878.137, 0.0, 0.0]  # Position on x-axis [km]
    v_init = [0.0, 7.546, 1.2]     # Velocity vector [km/s]

    initial_state = r_init + v_init

    # Propagate for roughly 1 orbital period (~5400 seconds)
    t_max = 5400.0
    t_eval = np.linspace(0, t_max, 500)

    # Run the propagator
    sol = propagate_orbit(initial_state, t_span=(0, t_max), t_eval=t_eval)

    # Display results
    print(f"\nPropagation Status: {sol.message}")
    print(f"Number of evaluation points: {len(sol.t)}")
    print(f"Initial Position [km]: {initial_state[0:3]}")
    print(f"Final Position [km]:   {sol.y[0:3, -1]}")
    print(f"Final Velocity [km/s]: {sol.y[3:6, -1]}")

if __name__ == "__main__":
    main()