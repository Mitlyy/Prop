import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.fftpack import *
import pandas as pd
from matplotlib.animation import ArtistAnimation

# Constants
I0 = 1E+11
lambda0 = 0.3 * 1E-3
tau = 5 * 1E-12
c = 3 * 1E+8



# Time
T = 80000 * 1E-15
N = 16384*2  # dots
dt = T/N
t = np.arange(-T/2, T/2, dt)  # array


# Signal params
f0 = c/lambda0  # Central freq
omega0 = 5.6 * 1e+12
df = 1/T
f = np.arange(-N*df/2, N*df/2, df)
omega = 2 * math.pi * f

# water
n2_cr = 8.572 * 1e-13
g_cr = omega0 * n2_cr / c

# ethanol
n_cr = pd.read_csv('n_eth.csv', delimiter = ",", header = None)[0].to_numpy()
k_z_cr = omega * n_cr / c

# Pulse params
E0 = math.sqrt(I0)
e_tau = np.exp(-2*(t/tau)**2)
E_in = E0 * (e_tau * np.sin(omega0 * t))

# Prop params
z_cr = 0.001  # crystal thickness

n2_cr1 = 5.05 * 1E-14
L_nl = lambda0/(16 * n2_cr1 * I0)
temp =  round(z_cr/(L_nl*0.01))
h_cr = z_cr/temp


freq = f
frames = []
frames2 = []
fig, ax = plt.subplot_mosaic([['a', 'a'], ['b', 'b']])
print(max((g_cr * (abs(E_in) ** 2) * h_cr)))
def split_step(f, steps):
	if steps != 0:
		D = 1  # fftshift(np.exp(-1j*k_z_cr*h_cr/2))
		N = np.exp(1j * g_cr * (abs(f) ** 2) * h_cr)
		y = fftshift(fft(f))
		y[y.shape[0]//2:] = (y[:y.shape[0]//2])[::-1]
		y = y * D
		y = ifft(y)
		y = y * N
		y = fftshift(fft(y))
		y[y.shape[0] // 2:] = (y[:y.shape[0] // 2])[::-1]
		y = y * D
		y = ifft(y)

		# Animation
		line, = ax['a'].plot(t, y, color = 'black')
		line2, = ax['b'].plot(freq, abs(fftshift(fft(y))), color = 'black')
		frames.append([line, line2])

		###
		return split_step(y, steps-1)
	else:
		return f

Eout = split_step(E_in, temp)



ax['a'].plot(t, E_in, dashes=[4, 2])
# ax['a'].plot(t, Eout, color="black")
ax['a'].legend(['In', 'Out'], loc="lower right")
ax['a'].set_xlim([-0.5 * 1E-11, 0.7 * 1E-11])

ax['b'].plot(f, abs(fftshift(fft(E_in))), dashes=[4, 2])
# ax['b'].plot(f, abs(fftshift(fft(Eout))), color = 'black')
ax['b'].legend(["In", "Out"], loc="upper right")
ax['b'].set_xlim([0, 4 * 1E+12])
# ax['b'].set_ylim([0, 0.25 * 1E+8])


animation = ArtistAnimation(fig, artists = frames, interval = 200, blit = True, repeat = False)
# animation2 = ArtistAnimation(fig, artists = frames2, interval = 200, blit = True, repeat = False)

animation.save('1.mp4', fps = 7, extra_args = ['-vcodec', 'libx264'])
plt.show()
