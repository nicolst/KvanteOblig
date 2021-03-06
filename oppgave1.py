import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

hbar = 1.055E-34 #Plancks reduserte konstant
m = 9.10938356E-31 #Elektronmassen
#k0 = 1.0E-30 # Hva i huleste
k0 = 1.0E10
dx = 1.0E-8 #
Ntot = 4000

d = hbar**2/(m*dx**2)
e = -hbar**2/(2*m*dx**2)
H = [[0]*Ntot for n in range(Ntot)]

for i in range(Ntot):
    print("Hamiltonian: {0} av {1}".format(i+1, Ntot))
    for j in range(Ntot):
        if i == j:
            H[i][j] = d
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


t_space = np.linspace(0, 1E-9, 200)

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

# plt.plot(t_space, analytical_uncertainty, '--')
# plt.plot(t_space, numerical_uncertainty, '-')
#
#
# for rho in rho_ts:
#     plt.plot(x, rho)
#
# plt.show()


fig = plt.figure('Wave packet animation', figsize=(16, 8))
ymax = np.max(rho_ts[0])
ax = plt.axes(xlim=(0, Ntot*dx), ylim=(0, ymax))
line, = ax.plot([], [], lw=1)






def init():
    line.set_data([], [])
    return line,


tidssteg = 1.0E-12


def animate(i):
    print(i)
    rho_t = rho_ts[i]
    #print(rho_t[len(rho_t)//2])
    line.set_data(x, rho_t)
    return line,

#plt.plot(x, [0]*len(x))
plt.xlabel(r"$x$ (m)", fontsize=20)
anim = animation.FuncAnimation(fig, animate, init_func=init, repeat=True, frames=len(t_space), interval=10, blit=True)



plt.show()


