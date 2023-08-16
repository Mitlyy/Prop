import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.fftpack import *
from matplotlib.animation import ArtistAnimation
import scienceplots

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])


### time
tau = 10e-12
t = np.linspace(-tau, tau, 10000)
f = np.linspace(-200*np.pi, 200*np.pi, 10000)

t_ = t / tau

# Input Pulse
I0 = 1e+11
E0 = math.sqrt(I0)
E_in = E0 * t_ * np.exp(-t_ ** 2)




mu_3 = 16 * 1e-9
mu2 = 2.11e+4/3e+8
temp = 40
h_cr = 0.005 / temp
# print(mu_3)
def split_step(f, steps, m):
	if steps != 0:
		y = fftshift(fft(f))
		y[y.shape[0] // 2:] = (y[:y.shape[0] // 2])[::-1]
		y = ifft(y)

		if m == "1":
			N1 = np.exp(1j * mu_3 * (abs(f)**2) * h_cr)
			print(N1)
		else:
			N1 = np.exp(1j * mu2 * (abs(f) ** 2) * h_cr)
		# N2 = np.exp(1j * g_cr2 * (abs(f) ** 2) * h_cr)
		# y = f
		y = y * N1
		# print(steps)

		return split_step(y, steps-1, m = m)
	else:
		return f

# E_out = split_step(E_in, temp, m = "2")
E_out2 = split_step(E_in, temp, m = "1")

Efft = abs(fftshift(fft(E_in)))
f_E_max = [f[i] for i in range(len(t)) if Efft[i] == max(Efft) and f[i] > 0]

plt.plot(f/f_E_max, abs(fftshift(fft(E_in))))
plt.plot(f/f_E_max, abs(fftshift(fft(E_out2))))
# plt.plot(t, E_in)
plt.xlim([0, 11])

plt.show()
