import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.fftpack import *
import pandas as pd
from matplotlib.animation import ArtistAnimation
import scienceplots
from PhysicalQuantities import q
import PhysicalQuantities.numpywrapper as nw

plt.style.use(['science', 'notebook', 'grid'])

### constants
dots = 10000
c = 3 * 1e+8          # m/s
l0 = 800e-9             # m
w_0 = 2 * np.pi * c / l0         # Hz
print(w_0)

n = 2.9e-16 #cm^2/W



### time
time_ = 80e-15
tau = 10e-15

t = np.linspace(-time_ / 2, time_ / 2, dots)
temp = 20
h_cr = 0.001 / temp


###
class material(object):
	def __init__(self, n2, chi2, N0):
		self.n2 = n2
		self.chi2 = chi2
		self.N0 = N0


class Pulse(object):
	def __init__(self, t, tau, n1 = 0, n2 = 0, N0 = 0, I = 0, chi2 = 180e-12):
		###

		self.time = t
		self.tau = tau
		self.n1 = n1
		self.n2 = n2
		self.mu3 = 2 * n2 * I / N0
		self.k1 = self.mu3 / c * w_0
		# self.k2 = chi2 * np.sqrt(I*10e+4)/N0/(3e+8 * tau)
		# print("K1",self.k1,'k2' , self.k2)
		self.dEdz = 0
		# print(self.k1)
		###
		self.A = None
		self.B = None
		self.nrm = None

		###
		self.E_out = None
		self.E0_k = None
		self.E_in = np.exp(-2 * (self.time/self.tau) ** 2) * np.sin(w_0 * self.time) #* q.V/q.m  # * math.sqrt(I)
		self.E0_t = np.copy(self.E_in)
		self.E0_test = np.copy(self.E_in)
	###
	def compute_fft(self):
		self.E0_k = fftshift(fft(self.E0_t))

	def compute_ifft(self):
		self.E0_t = ifft(self.E0_k)

		return None
	def plot(self):
		ax[1][0].plot(self.time, self.E0_test)
		self.E0_test = self.E0_test * (1 + self.E0_test**2 * h_cr * self.k1._ *\
		         ((-4*self.time/(self.tau**2)) + (np.cos(w_0 * self.time)/np.sin(w_0 * self.time))* w_0))
		ax[1][1].plot(self.time, self.E0_test/max( self.E0_test))
		self.E0_test = self.E0_test * (1 + self.E0_test ** 2 * h_cr * self.k1._ * \
		                               ((-4 * self.time / (self.tau ** 2)) + (
					                               np.cos(w_0 * self.time) / np.sin(w_0 * self.time)) * w_0))
		ax[2][0].plot(self.time, self.E0_test/max( self.E0_test))

	def compute_coeffs(self):
		# self.A = (1 + self.E0_t**2 * h_cr * self.k1._ *\
		#          ((-4*self.time/(self.tau**2)) + (np.cos(w_0 * self.time)/np.sin(w_0 * self.time))* w_0))
		self.A = np.exp(1j * h_cr * (self.k1._ * (abs(self.E0_t) ** 2)))
		self.B = None

	def norm(self):
		#self.nrm = [t[i] for i in range(len(t))[len(t) // 2:len(t) - 0*len(t) // 1] if
		           # abs(self.E0_k)[i] == max(abs(self.E0_k)) and t[i] > 0]

		self.nrm = t[abs(self.E0_k) == max(abs(self.E0_k))]
		return  self.nrm

	def split_step(self, steps):
		for i in range(steps):
			self.compute_fft()
			self.E0_k[self.E0_k.shape[0] // 2:] = (self.E0_k[:self.E0_k.shape[0] // 2])[::-1]
			self.compute_ifft()
			self.compute_coeffs()
			self.E0_t *= self.A
			self.E0_t /= max(self.E0_t)
			if i == steps - 1:
				self.compute_fft()
				self.E_out = self.E0_k
		else:
			pass
		return



ZnTe = material(n2 = 2.5e-12 * q.cm / q.W, chi2 = 180e-12 * q.pm / q.V, N0 = 1)


fig, ax = plt.subplots(3, 2, figsize = (8, 10))
f = Pulse(t, tau = tau, n2 = ZnTe.n2, N0 = ZnTe.N0, I = 3 * 10 ** 7 * q.W / q.cm**2)

f.split_step(temp)

f_max = max(f.norm())
# print(f_max)

def graphs():
	for i in range(5):
		start = 3
		if i < 3:
			f = Pulse(t, tau = tau, n2 = ZnTe.n2, N0 = ZnTe.N0, I = 3 * 10 ** (start + i)* q.W / q.cm**2)
			f.split_step(temp)
			ax[i][0].set_xlim([0, 6])
			ax[i][0].plot(t / f_max, abs(f.E_out))
			ax[i][0].legend([r"$3\cdot 10^{%s}$" % str(start + i)])
		else:
			f = Pulse(t, tau = tau, n2 = ZnTe.n2, N0 = ZnTe.N0, I = 3 * 10 ** (start + i)* q.W / q.cm**2)
			f.split_step(temp)
			ax[i - 3][1].set_xlim([0, 6])
			ax[i - 3][1].plot(t / f_max, abs(f.E_out))
			ax[i - 3][1].legend([r"$3\cdot 10^{%s}$" % str(start + i)])

graphs()
ax[2][1].plot(t, f.E_in)
# f.plot()
plt.show()
