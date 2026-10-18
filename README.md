# Spacecraft Mission Analysis & Flight Dynamics Portfolio

A modular, high-precision Python library for spacecraft trajectory propagation, orbital mechanics conversions, and maneuver design. Built to demonstrate core competencies in astrodynamics, numerical methods, and clean software architecture for mission analysis and flight dynamics engineering.

---

## Repository Structure

```text
mission-analysis-portfolio/
├── src/
│   ├── propagation/
│   │   ├── __init__.py
│   │   ├── propagator.py      # Numerical ODE integrator with J2 perturbations
│   │   ├── frames.py          # Classical Orbital Elements (COE) <-> ECI Cartesian conversions
│   │   ├── anomalies.py       # Kepler's equation solvers and anomaly conversions (Mean, Eccentric, True)
│   │   └── plot_orbit.py      # 3D orbit trajectory & Earth wireframe visualization
│   └── maneuvers/
│       ├── __init__.py
│       └── transfers.py       # Hohmann transfers & Universal Variable Lambert solver
├── README.md
└── requirements.txt
```

---

## Key Features & Modules

### 1. Core Propagation & Perturbations (`src/propagation/propagator.py`)
* Implements numerical orbit propagation using SciPy's high-order adaptive step-size IVP solvers (`solve_ivp`).
* **J2 Perturbation Modeling:** Accounts for Earth's oblate geometry (equatorial bulge), capturing critical real-world phenomena like **nodal regression** (essential for Sun-Synchronous Orbits) and **apsidal rotation**.

### 2. Coordinate Frames & Transformations (`src/propagation/frames.py`)
* Converts seamlessly between **Classical Orbital Elements (COEs)**—semi-major axis ($a$), eccentricity ($e$), inclination ($i$), RAAN ($\Omega$), argument of perigee ($\omega$), and true anomaly ($\nu$)—and Earth-Centered Inertial (ECI) Cartesian state vectors.

### 3. Orbit Anomalies & Kepler's Equation Solver (`src/propagation/anomalies.py`)
* Computes conversions between **Mean Anomaly ($M$)**, **Eccentric Anomaly ($E$)**, and **True Anomaly ($\nu$)**.
* Implements the **Newton-Raphson iterative method** to solve Kepler's transcendental equation for eccentric orbits with high precision.

### 4. Orbital Maneuvers & Boundary-Value Solver (`src/maneuvers/transfers.py`)
* **Impulsive Transfers:** Computes $\Delta V$ costs and time of flight for coplanar orbital changes (Hohmann transfers).
* **Lambert's Problem Solver:** A robust implementation of the Universal Variable formulation for solving boundary-value problems ($\mathbf{r}_1 \rightarrow \mathbf{r}_2$ given time-of-flight), crucial for rendezvous and interplanetary trajectory design.



## Tech Stack
* **Language:** Python 3.10+
* **Numerical Core:** NumPy, SciPy (`solve_ivp`, ODE integration, root-finding)
* **Visualization:** Matplotlib (3D wireframe projections)

