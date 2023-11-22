import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.fftpack import *
import pandas as pd
import scienceplots


plt.style.use(['science', 'notebook', 'grid'])


#Constants

dots = 2**13

# n = 2.9e-16                     #cm^2/W

time_ = 80e-12
# tau = 10e-12

# t = np.linspace(-tau_ / 2, time_ / 2, dots)      #ПОМЕПНЯТЬ ВРЕМЯ ИНТЕГРИРОВАНИЯ
steps = 20
h_cr = 0.005 / steps
c = 3 * 1e+8                    # m/s

  
class Pulse(object):
    def __init__(self, tau, n1 = 0, n2 = 0, N0 = 0, I = 0, l = 0):
     
        ###
        self.time = np.linspace(-tau *20, tau*20, dots)   
        self.tau = tau
        self.n_quad = n1
        self.n_cube = n2
        
        self.l0 = l                    # m
        self.w_0 = 2 * np.pi * c / self.l0        # Hz
        
        #cube
        self.mu2 = self.n_quad * np.sqrt(I) / N0
        self.k1 = self.mu2 / c * self.w_0
        
        #quad
        self.mu3 = 2 * self.n_cube * I / N0
        self.k2 = self.mu3 / c * self.w_0
        print("K1:", self.k1, "\n", "K2:", self.k2)
        # self.k2 = chi2 * np.sqrt(I*10e+4)/N0/(3e+8 * tau)
        # print("K1",self.k1,'k2' , self.k2)
  
  
        self.dEdz = 0
        ###
        self.A_split = None
        # self.B = None
        self.nrm = None

        ###
        self.E_out = None
        self.E0_k = None
        self.E_in = 0
        # self.E_in = np.exp(-2 * (self.time/self.tau) ** 2) * np.sin(w_0 * self.time) #* q.V/q.m  # * math.sqrt(I)
        
  
    ###
    def gauss(self, t_):
        ans = np.exp(-2 * (t_/self.tau) ** 2) * np.sin(self.w_0 * t_)
        return ans
    
        
    def compute_fft(self):
        self.E0_k = fftshift(fft(self.E0_t))

    def compute_ifft(self):
        self.E0_t = ifft(self.E0_k)

    def compute_coeffs(self):
        self.A_split = np.exp(1j * h_cr * (self.k1 * self.E0_t + self.k2 * self.E0_t**2))

    def norm(self):
        #self.nrm = [t[i] for i in range(len(t))[len(t) // 2:len(t) - 0*len(t) // 1] if
                   # abs(self.E0_k)[i] == max(abs(self.E0_k)) and t[i] > 0]

        self.nrm = self.time[:len(self.E0_k)][abs(self.E0_k) == max(abs(self.E0_k))]
        print("NOrm:", self.nrm)
        return  self.nrm

    def split_step(self, steps):
        self.E_in = self.gauss(self.time)
        self.E0_t= np.copy(self.E_in)
        print(self.k1, self.k2)
        
        for i in range(steps):
            self.compute_fft()
            self.E0_k[self.E0_k.shape[0] // 2:] = (self.E0_k[:self.E0_k.shape[0] // 2])[::-1]
            self.compute_ifft()
            self.compute_coeffs()
            self.E0_t *= self.A_split
            if i == steps - 1:
               
                self.compute_fft()
                self.E0_k /= max(self.E0_k)
                self.E_out = self.E0_k
        else:
            pass
        
    def crank_nik(self, steps):
        self.E_in = self.gauss(self.time)
        self.E0_t = self.crank(self.E_in, steps)
        

        self.compute_fft()
        self.E0_k /= max(self.E0_k)
        self.E_out = self.E0_k
    
    def crank(self, E_, steps):
        steps-=1
        # print(E_)
        E_h = E_[1:] + h_cr * np.diff(E_) * self.k2 * E_[1:] **2
        E_hk = 0
        # plt.plot(t[len(t)-len(E_h):], E_h)
        # plt.show() 
        # while not (abs(E_hk - E_h) < 1e-4).all():
        for i in range(10):
    
            # print(steps)
            E_h2 = 1/2 * (E_[len(E_)-len(E_h):] + E_h)
            E_hk = E_h[1:] + h_cr * np.diff(E_h2) * self.k2 * E_h[1:] **2
            
            E_h = E_hk
        if steps == 0:
            return E_hk
        return self.crank(E_hk, steps)


    # def crank_nik(self, steps):
    #     self.E_in = self.gauss(self.time)
    #     self.E0_t= np.copy(self.E_in)
    #     # for i in range(steps):
    #     # self.B1 = (1 + self.E0_t**2 * h_cr * self.k1 *\
    #     #     ((-4*self.time/(self.tau**2)) + (np.cos(self.w_0 * self.time)/np.sin(self.w_0 * self.time))* self.w_0))
        
    #     self.B2 = np.diff(self.E_in)


    #     plt.plot(self.time[1:]/self.tau, self.B2)

    #     # self.E0_t /= max(self.E0_t)
    #     # self.compute_fft()
    #     # self.E_out = self.E0_k

        
        
            
            

class material(object):
    def __init__(self, n2, n1, N0):
        self.n2 = n2
        self.n1 = n1
        self.N0 = N0
  