# Lunar Transfer Orbit (LTO) Mission Analysis Module

A specialized astrodynamics module designed to compute cislunar trajectories, Trans-Injection Delta-V requirements, and comparative Earth-Moon transfer mechanics. This module extends the baseline portfolio from Low Earth Orbit (LEO) operations into deep-space and interplanetary mission design.

---

## Directory Structure

```text
lunar_transfer/
├── README.md               # Module specifications and engineering overview
├── lto_calculator.py       # TLI delta-V and transfer orbit propagator 
└── test_lto.py             # Validation scripts and boundary-value verification 
```
## Mission Specifications & Parameters
1. Initial State & Parking Orbit
   Orbit Type: Circular Low Earth Orbit (LEO) parking orbit.
   Altitude ($h_{\text{park}}$): 300.0 km 
   Equatorial Radius of Earth ($R_\oplus$): 6378.137 km
   Parking Radius ($r_{\text{park}}$): 6678.137 km
   Earth Gravitational Parameter ($\mu_\oplus$): $398600.4418 \text{ km}^3/\text{s}^2$2. 
2. Target Body (The Moon)
   Semi-Major Axis ($a_{\text{moon}}$): 384,400.0 km
   Lunar Gravitational Parameter ($\mu_{\text{moon}}$): $4902.8 \text{ km}^3/\text{s}^2$

## Core Mathematical Formulations
   Parking Orbit Velocity
   The baseline velocity of the spacecraft in its circular parking orbit prior to injection:
   $$v_{\text{park}} = \sqrt{\frac{\mu_\oplus}{r_{\text{park}}}}$$
   Transfer Orbit 
   Semi-Major AxisThe transfer ellipse is engineered such that perigee aligns with the LEO parking radius and apogee reaches the orbital radius of the Moon:
   $$a_{\text{trans}} = \frac{r_{\text{park}} + r_{\text{moon}}}{2}$$
   Perigee Velocity and Trans-Lunar Injection (TLI)
   Delta-VUsing the vis-viva energy conservation equation evaluated at perigee, the velocity required on the transfer ellipse ($v_{\text{perigee}}$) determines the impulsive burn magnitude ($\Delta v_{\text{TLI}}$):
   $$v_{\text{perigee}} = \sqrt{\mu_\oplus \left( \frac{2}{r_{\text{park}}} - \frac{1}{a_{\text{trans}}} \right)}$$
   $$\Delta v_{\text{TLI}} = v_{\text{perigee}} - v_{\text{park}}$$
   Transfer Time of FlightBased on Kepler's Third Law for half of the elliptical transfer orbit:
   $$t_{\text{tof}} = \pi \sqrt{\frac{a_{\text{trans}}^3}{\mu_\oplus}}$$
   