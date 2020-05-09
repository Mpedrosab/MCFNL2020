import numpy as np
import math
import matplotlib.pyplot as plt


class MeasureTransmittance:
    '''
        Measure transmittance through a layer
    '''
    def __init__(self,layer,probe):
        self._initIndex = layer.indices[0]
        self._endIndex = layer.indices[-1]

        #Get time
        self.t = np.array(probe['time'])

        #Get E field at both sides of the layer
        self._initE = np.zeros(self.t.size)
        self._endE = np.zeros(self.t.size)
        i=0
        for p in probe['values']:

            if i==0:        #First element is a list with an array inside (do not know why)
                self._initE[i] = p[0][self._initIndex ]
                self._endE[i] = p[0][self._endIndex ]
            else:
                self._initE[i] = p[self._initIndex ]
                self._endE[i] = p[self._endIndex ]               
            i+=1

    def getT(self, initE, endE):
        self.T = endE-initE
        return self.T

    def AmplVsFreq(self):
        self.T = self.getT(self._initE, self._endE)
        self.f_Tq = np.fft.fftfreq(len(self.t)) / (self.t[1]-self.t[0])
        self.Tq = np.fft.fft(self.T )

        return (self.f_Tq,np.abs(self.Tq))

    def PlotTransmittance(self,freq,Trans):
        plt.figure()
        plt.plot(freq, np.abs(Trans))
        plt.figure()
        plt.plot(np.fft.fftshift(freq), np.fft.fftshift(np.abs(Trans)))
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Transmittance')
        plt.show()
