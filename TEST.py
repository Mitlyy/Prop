import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import *

c = 3e+8
l0 = 900e-9
dots = 2**12
time_ = 80e-12
tau = 10e-16
w_0 = 2 * np.pi * c / l0
h_cr = 0.005 / 20

t = np.linspace(-tau * 10, tau *10, dots)

def gauss(t_):
    ans = np.exp(-2 * (t/tau) ** 2) * np.sin(w_0 * t_)
    return ans
plt.plot(t, gauss(t))
plt.show()

def crank(E_, steps):
    steps-=1
    print(E_)
    E_h = E_[1:] + h_cr * np.diff(E_) * 10
    E_hk = 0
    plt.plot(t[len(t)-len(E_h):], E_h)
    plt.show() 
    # while not (abs(E_hk - E_h) < 1e-4).all():
    for i in range(10):
 
        print(steps)
        E_h2 = 1/2 * (E_[len(E_)-len(E_h):] + E_h)
        E_hk = E_h[1:] + h_cr * np.diff(E_h2) * 10
        
        E_h = E_hk
    if steps == 0:
        return E_hk
    return crank(E_hk, steps)
E_in = gauss(t)
E0_t = crank(E_in, 20)


E_out = fftshift(fft(E0_t))
plt.plot(t[len(t)-len(E_out):], abs(E_out))
plt.show()
    