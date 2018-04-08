import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

hbar = 1.055E-34

class Partikkel:
    def __init__(self, m, k0, dx, V0, N, sigma, x_offset=0):
        self.m = m
        self.k0 = k0
        self.dx = dx
        self.V0 = V0
        self.N = N
        self.sigma = sigma

        V = [V0] * 4 * N + [V0 * ((n - N) / (N * 1.0)) ** 2 for n in range(2 * N + 1)] + [V0] * 4 * N
        V = np.asarray(V)
        self.V = V

        self.Ntot = len(V)

        self.d = np.array([v + hbar ** 2 / (m * dx ** 2) for v in V])
        self.e = -hbar ** 2 / (2 * m * dx ** 2)

        energy, psi_matrix = la.eigh_tridiagonal(self.d, np.array([self.e] * (self.Ntot - 1)))
        self.energy = energy
        self.psi_matrix = psi_matrix
        self.psi_matrix_complex = psi_matrix*(1.0 + 0.0j)

        self.x = np.asarray([dx * n for n in range(self.Ntot)])
        self.x0 = self.x[len(self.x) // 2] + x_offset

        normfactor = (1.0 + 0.0j) * (2 * np.pi * sigma ** 2) ** (-0.25)
        gaussinit = (1.0 + 0.0j) * np.exp(-(self.x - self.x0) ** 2 / (4 * sigma ** 2))
        planewavefactor = np.exp(1j * k0 * self.x)
        self.Psi0 = normfactor * gaussinit * planewavefactor

        self.c = np.zeros(self.Ntot, dtype=np.complex128)
        for n in range(self.Ntot):
            print("Dotting: {0} out of {1}".format(n + 1, self.Ntot))
            self.c[n] = np.vdot(self.psi_matrix_complex[:, n], self.Psi0)

    def calculate(self, max_t, steps):
        self.t_space = np.linspace(0, max_t, steps)

        # Analytisk usikkerhet
        self.analytical_uncertainty = np.sqrt(self.sigma ** 2 + hbar ** 2 * self.t_space ** 2 / (4 * self.m ** 2 * self.sigma ** 2))

        # Numerisk usikkerhet
        self.numerical_uncertainty = []
        self.rho_ts = []
        for i in range(len(self.t_space)):
            t = self.t_space[i]
            print("T: {0} av {1}".format(i+1, len(self.t_space)))
            Psi_t = np.zeros(self.Ntot, dtype=np.complex128)
            for n in range(self.Ntot):
                Psi_t = Psi_t + self.c[n] * self.psi_matrix_complex[:, n] * np.exp(-1j * self.energy[n] * t / hbar)

            rho_t = np.abs(Psi_t) ** 2
            self.rho_ts.append(rho_t)
            x2mean = self.dx * np.dot(self.x ** 2, rho_t)
            xmean2 = (self.dx * np.dot(self.x, rho_t)) ** 2
            deltax = np.sqrt(x2mean - xmean2)

            self.numerical_uncertainty.append(deltax)






p1 = Partikkel(9.1E-31, 1.0E10, 1.0E-10, 0.0, 100, 1.0E-9)
print(p1.psi_matrix)



























