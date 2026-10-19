
import numpy as np
from propagation.anomalies import Eccentric_anomaly, true_anomaly

M = 1.0
e = 0.1
tol = 1e-10
max_iter = 100

E = Eccentric_anomaly(M, e, tol, max_iter)
theta = true_anomaly(E, e)

print("\n--- Results ---")
print(f"Mean anomaly      : {M:.10f} rad")
print(f"Eccentric anomaly : {E:.10f} rad")
print(f"True anomaly      : {theta:.10f} rad")
print(f"True anomaly      : {np.degrees(theta):.6f} degrees")

# Independent verification of Kepler's equation
residual = E - e * np.sin(E) - M
print(f"Final residual    : {residual:.3e}")

assert abs(residual) <= tol, "Kepler equation verification failed!"
print("TEST PASSED")
