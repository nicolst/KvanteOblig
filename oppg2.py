import partikkel
import matplotlib.pyplot as plt
import numpy as np

m_e = 9.11E-31
hbar = 1.055E-34
k0 = 1.0E8
V0 = 3*hbar**2 * k0**2 / m_e
dx = 1.0E-10
sigma = 100*dx
time_step = 1.0E-13

p1 = partikkel.Partikkel(m_e, k0, dx, V0, 400, sigma, x_offset=0)
p1.calculate(time_step, 900)

plt.figure(1)
plt.title(r"Usikkerhet i posisjon")
plt.xlabel(r"$x$ / (m)", fontsize=20)
plt.ylabel(r"$\Delta x$ / m", rotation=0, fontsize=20)
plt.plot(p1.t_space, p1.numerical_uncertainty, 'k-', label="Usikkerhet")
#plt.plot(p1.t_space, p1.analytical_uncertainty, 'k--', label="Analytisk")
plt.legend()
plt.show()

# Usikkerhet i impuls
Psi_t_behind = hbar/1j * p1.Psi_ts
Psi_t_ahead = np.asarray([0] + list(Psi_t_behind[:-1]))
Psi_t_diff = (Psi_t_ahead - Psi_t_behind)/time_step

Psi_t_behind_2 = hbar/1j * Psi_t_diff
Psi_t_ahead_2 = np.asarray([0] + list(Psi_t_behind_2[:-1]))
Psi_t_diff_2 = (Psi_t_ahead_2 - Psi_t_behind_2)/time_step

rho_impulser = np.vdot(p1.psi_matrix_complex, Psi_t_diff)
rho_impulser_2 = np.vdot(p1.psi_matrix_complex, Psi_t_diff_2)

p2mean_m = dx * rho_impulser_2
pmean2_m = (dx * rho_impulser)**2
deltap = p2mean_m - pmean2_m

plt.figure(2)
plt.title(r"Usikkerhet i impuls")
plt.xlabel(r"$p$ / (kg m / s)", fontsize=20)
plt.ylabel(r"$\Delta p$ / (kg m / s)", rotation=0, fontsize=20)
plt.plot(p1.t_space, p1.numerical_uncertainty, 'k-', label="Usikkerhet")
#plt.plot(p1.t_space, p1.analytical_uncertainty, 'k--', label="Analytisk")
plt.legend()
plt.show()









