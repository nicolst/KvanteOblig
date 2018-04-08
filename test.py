from matplotlib import animation
import partikkel
import matplotlib.pyplot as plt
import numpy as np

m_e = 9.1E-31
hbar = 1.055E-34
k0 = 1.0E8
V0 = 3*hbar**2 * k0**2 / m_e
dx = 1.0E-10
sigma = 100*dx

p1 = partikkel.Partikkel(m_e, k0, dx, V0, 400, sigma, x_offset=0)
p1.calculate(1.0E-13, 900)

fig = plt.figure('Wave packet animation', figsize=(16, 8))
ymax = np.max(p1.rho_ts[0])
ax = plt.axes(xlim=(0, p1.Ntot*p1.dx), ylim=(0, ymax))
line, = ax.plot([], [], lw=1)
ax.plot(p1.x, p1.V * (ymax / 2 / p1.V0))

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
plt.title("")
plt.plot(p1.t_space, p1.analytical_uncertainty, 'k--')
plt.plot(p1.t_space, p1.numerical_uncertainty, 'k-')

#writer = animation.FFMpegWriter(fps=60, bitrate=3000)
#anim.save('test.mp4', writer=writer)

plt.show()