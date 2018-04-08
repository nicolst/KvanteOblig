from matplotlib import animation
import partikkel
import matplotlib.pyplot as plt
import numpy as np

m_e = 9.1E-31
hbar = 1.055E-34
k0 = 1.0E10
V0 = 3*hbar**2 * k0**2 / m_e

p1 = partikkel.Partikkel(9.1E-31, 1.0E10, 1.0E-7, V0, 400, 1.0E-6, x_offset=-5.0E-6)
p1.calculate(1.0E-7, 500)

fig = plt.figure('Wave packet animation', figsize=(16, 8))
ymax = np.max(p1.rho_ts[0])
ax = plt.axes(xlim=(0, p1.Ntot*p1.dx), ylim=(0, ymax))
line, = ax.plot([], [], lw=1)


def init():
    line.set_data([], [])
    return line,


def animate(i):
    print(i)
    rho_t = p1.rho_ts[i]
    line.set_data(p1.x, rho_t)
    fig.suptitle(p1.numerical_uncertainty[i])
    return line,

plt.xlabel(r"$x$ (m)", fontsize=20)
anim = animation.FuncAnimation(fig, animate, init_func=init, repeat=True, frames=len(p1.t_space), interval=17, blit=True)

plt.figure(2)


plt.show()