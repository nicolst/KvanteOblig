import partikkel
import matplotlib.pyplot as plt

m_e = 9.11E-31
hbar = 1.055E-34
k0 = 1.0E8
V0 = 3*hbar**2 * k0**2 / m_e
dx = 1.0E-10
sigma = 100*dx

p1 = partikkel.Partikkel(m_e, k0, dx, 0, 400, sigma, x_offset=0)
p1.calculate(1.0E-13, 900)

plt.figure(1)
plt.title(r"Usikkerhet i posisjon")
