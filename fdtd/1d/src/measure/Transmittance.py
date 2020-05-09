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

            self._initE[i] = p[0][self._initIndex ]
            self._finalE = p[0][self._endIndex ]
            i+=1


        #Transmittance in time
        self.T = self.getT(self._initE,self._endE)

        #Transmittance in freq
        self.TFreq = self.AmplVsFreq(self.t,self.T)
    def getT(self, initE, endE):
        return endE-initE

    def AmplVsFreq(self,t,f):

        fq = np.fft.fftfreq(len(t)) / (t[1]-t[0])
        f_fq = np.fft.fft(f)

    '''
    plt.figure()
    plt.plot(t, f)

    plt.figure()
    plt.plot(np.fft.fftshift(fq), np.fft.fftshift(np.abs(f_fq)))

    plt.show()

    print('=== Program finished ===')
    '''