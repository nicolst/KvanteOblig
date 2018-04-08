import partikkel
import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fftshift

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
plt.ylabel(r"$\Delta x$ / m", fontsize=20)
plt.plot(p1.t_space, p1.numerical_uncertainty, 'k-', label="Usikkerhet")
#plt.plot(p1.t_space, p1.analytical_uncertainty, 'k--', label="Analytisk")
plt.legend()

Psi_ks = np.asarray([fftshift(psi) for psi in p1.Psi_ts])
Psi_ps = Psi_ks * hbar
rho_p = np.abs(Psi_ps)**2

deltaps = []

for i in range(900):
    p2mean_m = rho_p*k0
    pmean2_m = (rho_p*k0)**2
    deltap = p2mean_m - pmean2_m
    deltaps.append(deltap)

plt.figure(2)
plt.title(r"Usikkerhet i impuls")
plt.xlabel(r"$p$ / (kg m / s)", fontsize=20)
plt.ylabel(r"$\Delta p$ / (kg m / s)", fontsize=20)
plt.plot(p1.t_space, deltaps, 'k-', label="Usikkerhet")
#plt.plot(p1.t_space, p1.analytical_uncertainty, 'k--', label="Analytisk")
plt.legend()
plt.show()









