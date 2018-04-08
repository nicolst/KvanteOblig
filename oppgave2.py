import matplotlib.pyplot as plt
import numpy as np

hbar = 1.055E-34 #Plancks reduserte konstant
m = 9.10938356E-31 #Elektronmassen
#k0 = 1.0E-30 # Hva i huleste
k0 = 1.0E10
dx = 1.0E-8 #
N = 200
#V0=5 * k0**2 * hbar**2 / (2 * m)
V0=1.6E-19

V = [V0]*4*N + [V0*((n-N)/(N*1.0))**2 for n in range(2*N+1)] + [V0]*4*N
V = np.asarray(V)

Ntot = len(V)



d = [(v + hbar**2/(m*dx**2)) for v in V]
e = -hbar**2/(2*m*dx**2)
H = [[0]*Ntot for n in range(Ntot)]

for i in range(Ntot):
    print("Hamiltonian: {0} av {1}".format(i+1, Ntot))
    for j in range(Ntot):
        if i == j:
            H[i][j] = d[i]
        elif abs(i-j) == 1:
            H[i][j] = e

energy, psi_matrix = np.linalg.eigh(H) # Egenverdier (energier) og egenvektorer (diskret bølgefunksjon)

print(np.shape(psi_matrix), np.shape(H))

x = np.asarray([dx*n for n in range(Ntot)])

# Starttilstand
x0 = x[len(x)//2]
sigma = 100*dx
normfactor = (2*np.pi*sigma**2)**(-0.25)
gaussinit = np.exp(-(x-x0)**2/(4*sigma**2))
planewavefactor = np.exp(1j*k0*x)
Psi0 = normfactor*gaussinit*planewavefactor
print(Psi0[len(Psi0)//2])

# Utviklingskoeffisienter
psi_matrix_complex = psi_matrix*(1.0 + 0.0j)
c = np.zeros(Ntot, dtype = np.complex128)
for n in range(Ntot):
    print("Dotting: {0} out of {1}".format(n+1, Ntot))
    c[n] = np.vdot(psi_matrix_complex[:,n], Psi0)


t_space = np.linspace(0, 1E-12, 200)

# Analytisk usikkerhet
analytical_uncertainty = np.sqrt(sigma**2 + hbar**2 * t_space**2 / (4 * m**2 * sigma**2))

# Numerisk usikkerhet
numerical_uncertainty = []
rho_ts = []
for t in t_space:
    print("T: {0} av {1}".format(t, t_space[-1]))
    Psi_t = np.zeros(Ntot, dtype=np.complex128)
    for n in range(Ntot):
        Psi_t = Psi_t + c[n]*psi_matrix_complex[:, n]*np.exp(-1j*energy[n]*t/hbar)

    rho_t = np.abs(Psi_t)**2
    rho_ts.append(rho_t)
    x2mean = dx*np.dot(x**2, rho_t)
    xmean2 = (dx*np.dot(x, rho_t))**2
    deltax = np.sqrt(x2mean - xmean2)

    numerical_uncertainty.append(deltax)

plt.figure(1)
plt.plot(t_space, analytical_uncertainty, '--')
plt.plot(t_space, numerical_uncertainty, '-')

plt.figure(2)
plt.plot(x, V, 'k-')
plt.plot(x, rho_ts[0])

plt.figure(3)
plt.plot(x, V, 'k-')
plt.plot(x, rho_ts[len(rho_ts)//2])


plt.show()


