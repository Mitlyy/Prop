import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.fftpack import *
import pandas as pd
# Constants
I0 = 1E+11
lambda0 = 0.3 * 1E-3
tau = 2000 * 1E-15
c = 3 * 1E+8
eps = 8.85 * 1E-12


# Time
T = 80000 * 1E-15
N = 16384*2  # dots
dt = T/N
t = np.arange(-T/2, T/2, dt)  # array


# Signal params
f0 = c/lambda0  # Central freq
omega0 = 2 * math.pi * f0
df = 1/T
f = np.arange(-(N-1)*df/2, N*df/2, df)[:32768]
omega = 2 * math.pi * f

# Space params #
N0_air = 1
a_air = 2.0959 * 1E-45
b_air = 0

# water
n2_cr = 8.572 * 1e-13
g_cr = omega0 * n2_cr / c
# ethanol
n_cr = pd.read_csv('n_eth.csv', delimiter = ",", header = None)[0].to_numpy()
k_z_cr = omega0 * n_cr / c

# Pulse params
E0 = math.sqrt(I0)
e_tau = math.e**(-2*(t/tau)**2)
E_in = E0 * (e_tau * np.sin(omega0 * t))

# Prop params
F = 0.5  # focal length
z_dom = 2 * F
N_z = 20  # steps
dz = z_dom/N_z
z_cr = 0.001  # crystal thickness

N0_water = 1.8
n2_cr1 = 5.05 * 1E-14
L_nl = lambda0/(16 * n2_cr1 * I0)
temp = round(z_cr/(L_nl*0.01))
h_cr = z_cr/temp

D = fftshift(np.exp(-1j*k_z_cr[:N]*(h_cr/2)))
Ein = E_in.copy()

for i in range(temp):
	N = 1j * g_cr * abs(Ein)**2
	y = fft(Ein)
	y = D * y
	y = ifft(y)
	y = np.exp(h_cr * N) * y
	y = fft(y)
	y = D * y
	y = ifft(y)


fig, ax = plt.subplot_mosaic([['a', 'b'], ['a', 'b']])
ax['a'].plot(t, y.real)
ax['a'].plot(t, Ein, color="black")
ax['a'].legend(['In', 'Out'], loc="lower right")
ax['a'].set_xlim([-0.5 * 1E-11, 0.5 * 1E-11])

ax['b'].plot(f, abs(fftshift(fft(Ein))))
ax['b'].plot(f, abs(fftshift(fft(y))), dashes=[4, 2], color = 'black')
ax['b'].legend(["In", "Out"], loc="upper right")
ax['b'].set_xlim([0, 4 * 1E+12])
plt.show()
