import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light, epsilon_0, mu_0

class MeasureTransmittance:
    '''
        Measure transmittance through a layer
    '''
    def __init__(self,layer,t,field, fieldDispersed):

        #Border of the layer
        self._endIndex = layer.indices[-1]

        #Get time
        self.t = np.array(t)

        #Get E field at both sides of the layer
        self._initE = np.zeros(self.t.size)
        self._endE = np.zeros(self.t.size)

        for i in range(0,len(field),1):

            if i==0:        #First element is a list with an array inside (do not know why)
                self._initE[i] = field[i][0][self._endIndex ]
                self._endE[i] = fieldDispersed[i][0][self._endIndex ]
            else:
                self._initE[i] = field[i][self._endIndex]
                self._endE[i] = fieldDispersed[i][self._endIndex ]               

        print(i)

    def getT(self, initE, endE):
        self.ratio = np.abs(np.fft.fft(endE) )/np.abs(np.fft.fft(initE ))
        #self.T = np.nan_to_num(self.T)
        return self.ratio

    def AmplVsFreq(self):
        self.Tq = self.getT(self._initE, self._endE)
        self.f_Tq = np.fft.fftfreq(len(self.t)) / (self.t[1]-self.t[0])


        return (self.f_Tq,np.abs(self.Tq))


class AnalyticTransmittance:

    def __init__(self, layer):
        self.thickness = layer.width
        self.epsilon_r = layer.epsilon
        self.ap=layer.ap
        self.cp=layer.cp
        self._eta_0 = np.sqrt(mu_0/epsilon_0)

        try:
            self.mu_r = layer.mu
        except:
            self.mu_r = 1.0
            Warning('mu not defined. Setting to 1.0')

        try:
            self.sigma = sigma
        except:
            self.sigma = 0.0
            Warning('sigma not defined. Setting to 0.0')

    def epsilon_c(self, omega):
        complexEps = 0
        for i in range(0,len(self.cp),1):
            complexEps += (self.cp[i]/(complex(0,1)*omega-self.ap[i])) + (np.conj(self.cp[i])/(complex(0,1)*omega-np.conj(self.ap[i]))) 
        return self.epsilon_r*epsilon_0 + epsilon_0 * complexEps

    def mu_c(self, omega):
        return self.mu_r * mu_0

    def gamma(self, omega):
        return complex(0,1) * omega * \
            np.sqrt(self.epsilon_c(omega) * self.mu_c(omega))

    def eta(self, omega):
        return np.sqrt(self.mu_c(omega) / self.epsilon_c(omega))

    def phi(self, omega):
        gd  = self.gamma(omega) * self.thickness
        eta = self.eta(omega)
        return np.array([[np.cosh(gd),      np.sinh(gd) * eta], \
                         [np.sinh(gd) /eta, np.cosh(gd)      ]])

    def _den(self, omega):
        phi = self.phi(omega)
        return phi[0,0]*self._eta_0 + phi[0,1] + phi[1,0]*self._eta_0**2 + phi[1,1]*self._eta_0
        
    def T(self, omega):
        sol = 2*self._eta_0 / self._den(omega)
        return  sol

    def R(self, omega): 
        phi = self.phi(omega)
        return \
            (phi[0,0]*self._eta_0 + phi[0,1] - phi[1,0]*self._eta_0**2 - phi[1,1]*self._eta_0) / \
            self._den(omega)

class PlotTransmittance:
    def __init__(self,freq,Trans):
        #plt.figure()
        #plt.plot(freq, np.abs(Trans))
        plt.figure()
        Trans = Trans[freq>=0]
        freq = freq[freq>=0]
        plt.plot(freq, Trans)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Transmittance')
        plt.show()
