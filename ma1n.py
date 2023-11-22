import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.fftpack import *
import pandas as pd
from matplotlib.animation import ArtistAnimation
import scienceplots

plt.style.use(['science', 'notebook', 'grid', 'Solarize_Light2'])

### constants
dots = 10000
c = 3 * 1e+8  # m/s
l0 = 800e-9 #m
w_0 = 2 * np.pi * c / l0  # s^-1


### time
tau = 1e-14
t = np.linspace(-tau * 20, tau * 20, dots)
temp = 20
h_cr = 0.005 / temp


###
class material(object):
	def __init__(self, n2, chi2, N0):
		self.n2 = n2
		self.chi2 = chi2
		self.N0 = N0


class Pulse(object):
	def __init__(self, t, tau, n1 = 0, n2 = 0, N0 = 0, I = 0, chi2 = 180e-12):
		###
		self.time = t / tau  # ms
		self.n1 = n1
		self.n2 = n2
		self.mu3 = 2 * n2 * I / N0
		self.k1 = self.mu3 / c * w_0
		self.k2 = chi2 * np.sqrt(I*10e+4)/N0/(3e+8 * tau)
		print("K1",self.k1,'k2' , self.k2)

		# print(self.k1)
		###
		self.A = None
		self.B = None
		self.nrm = None

		###
		self.E_out = None
		self.E0_k = None
		self.E_in = np.exp(-2 * self.time ** 2) * np.sin(w_0 * t)  # * math.sqrt(I)
		self.E0_t = np.copy(self.E_in)

	###

	def compute_fft(self):
		self.E0_k = fftshift(fft(self.E0_t))

	def compute_ifft(self):
		self.E0_t = ifft(self.E0_k)

	def compute_coeffs(self):
		self.A = np.exp(1j * h_cr * (self.k1 * (self.E0_t) ** 2)) #+ self.k2 * (abs(self.E0_t))) )
		self.B = None

	def norm(self):
		#self.nrm = [t[i] for i in range(len(t))[len(t) // 2:len(t) - 0*len(t) // 1] if
		           # abs(self.E0_k)[i] == max(abs(self.E0_k)) and t[i] > 0]
		self.nrm = t[abs(self.E0_k) == max(abs(self.E0_k))]
		return self.nrm

	def split_step(self, steps):
		for i in range(steps):
			self.compute_fft()
			self.E0_k[self.E0_k.shape[0] // 2:] = (self.E0_k[:self.E0_k.shape[0] // 2])[::-1]
			self.compute_ifft()
			self.compute_coeffs()
			self.E0_t *= self.A
			if i == steps - 1:
				self.compute_fft()
				self.E_out = self.E0_k
		else:
			pass
		return


ZnTe = material(n2 = 1.25e-12, chi2 = 180e-12, N0 = 2.79)
fig, ax = plt.subplots(3, 2, figsize = (8, 10))
f = Pulse(t, tau = tau, n2 = ZnTe.n2, N0 = ZnTe.N0, I = 3 * 10 ** 7)
f.split_step(temp)
f_max =abs(f.norm()[0])
print(f_max)
def graphs():
	for i in range(5):
		start = 3
		if i < 3:
			f = Pulse(t, tau = tau, n2 = ZnTe.n2, N0 = ZnTe.N0, I = 3 * 10 ** (start + i))
			f.split_step(temp)
			ax[i][0].set_xlim([0, 6])
			ax[i][0].plot(t / f_max, abs(f.E_out))
			ax[i][0].legend([r"$3\cdot 10^{%s}$" % str(start + i)])
		else:
			f = Pulse(t, tau = tau, n2 = ZnTe.n2, N0 = ZnTe.N0, I = 3 * 10 ** (start + i))
			f.split_step(temp)
			ax[i - 3][1].set_xlim([0, 6])
			ax[i - 3][1].plot(t / f_max, abs(f.E_out))
			ax[i - 3][1].legend([r"$3\cdot 10^{%s}$" % str(start + i)])

graphs()
ax[2][1].plot(t, f.E_in)

plt.show()