import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

hbar = 1.055E-34
m = 9.10938356E-31
k0 = 1.0E-10
dx = 1.0E-12
N = 200
Ntot = 10*N + 1

d = [hbar**2/(m*dx**2)]*Ntot
e = -hbar**2/(2*m*dx**2)
H = [[0]*Ntot for n in range(Ntot)]

for i in range(Ntot):
    print("Hamiltonian {0} of {1}".format(i, Ntot-1))
    for j in range(Ntot):
        if i == j:
            H[i][j] = d[i]
        elif abs(i-j) == 1:
            H[i][j] = e

H = np.asarray(H)
energy, psi_matrix = np.linalg.eigh(H)
print(np.shape(psi_matrix))

x = np.asarray([dx*n for n in range(Ntot)])
x0 = 5*N*dx
sigma = 300*dx
normfactor = (2*np.pi*sigma**2)**(-0.25)
gaussinit = np.exp(-(x-x0)**2/(4*sigma**2))
planewavefactor = np.exp(1j*k0*x)
Psi0 = normfactor*gaussinit*planewavefactor
#print(Psi0)

psi_matrix_complex = psi_matrix*(1.0 + 0.0j)
c = np.zeros(Ntot, dtype = np.complex128)
for n in range(Ntot):
    print("Dotting: {0} out of {1}".format(n, Ntot-1))
    c[n] = np.vdot(psi_matrix_complex[:,n], Psi0)

time_step = 3.0E-15
steps = 10
t_space = np.linspace(0, 1, 200)

analytical_uncertainty = np.sqrt(sigma**2 + (hbar*t_space/(2*m*sigma))**2)

time_psi = []
numerical_uncertainty = []
num = 1
for t in t_space:
    print("T_space: {0} out of {1}".format(num, len(t_space)))
    num += 1
    psi = np.dot(psi_matrix, c*np.exp(-1j*energy*t/hbar))
    rho_t = np.abs(psi)**2

    x2mean = dx*np.dot(x**2, rho_t)
    xmean2 = (dx*np.dot(x, rho_t))**2

    unc = np.sqrt(x2mean-xmean2)
    numerical_uncertainty.append(np.sqrt(x2mean-xmean2))


    time_psi.append(psi)


plt.plot(t_space, analytical_uncertainty, '-')
plt.plot(t_space, numerical_uncertainty, '--')

plt.show()


# for tim in time_psi:
#     plt.plot(x, tim)
#
# plt.show()


# fig = plt.figure('Wave packet animation', figsize=(16, 8))
# ymax = 1.0E3
# ax = plt.axes(xlim=(0, Ntot*dx), ylim=(0, ymax))
# line, = ax.plot([], [], lw=1)
#
#
#
#
#
#
# def init():
#     line.set_data([], [])
#     return line,
#
#
# tidssteg = 3.0E-20
#
#
# def animate(i):
#     print(i)
#     tt = i*tidssteg
#     Psi_t = np.zeros(Ntot, dtype=np.complex128)
#     for n in range(Ntot):
#         Psi_t += c[n]*psi_matrix_complex[:,n]*np.exp(-1j*energy[n]*tt/hbar)
#
#     rho_t = np.abs(Psi_t)**2
#     print(rho_t[len(rho_t)//2])
#     line.set_data(x, rho_t)
#     return line,
#
# plt.plot(x, [0]*len(x))
# plt.xlabel(r"$x$ (m)", fontsize=20)
# anim = animation.FuncAnimation(fig, animate, init_func=init, repeat=False, frames=2000, interval=500, blit=True)
#
#
# plt.show()














