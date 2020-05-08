import numpy as np
import copy
import math
import scipy.constants as sp

class DispersiveMedia:
    def __init__(self,mesh,layer):
        #Copy outside since it is needed in every loop. Speed up the code
        self._media=copy.deepcopy(layer)
        self._ap=np.empty(len(self._media["ap"]),dtype=complex)
        self._cp=np.empty(len(self._media["cp"]),dtype=complex)
        self._epsilon=self._media['permittivity'] 
        self.pos = self._media['startPosition']
        self.width = self._media['width']
        for i in range(0,len(self._media["ap"])):
            self._ap[i]=complex(self._media["ap"][i])
            self._cp[i]=complex(self._media["cp"][i])
        
        if (self._media['unitsFreq']!="Hz"):
            if (self._media['unitsFreq']=="eV"):
                self._ap = self._ap * sp.e/sp.h
                self._cp = self._cp * sp.e/sp.h
            else:
                raise ValueError("Invalid frequency units. Frequency must be in Hz or eV")
        #self._cp=self._cp * sp.epsilon_0
        # self._epsilon= self._epsilon * sp.epsilon_0
    #WATCH OUT! Normalize permittivity0
        #self._ap = self._ap / self._media['permittivity']
        #self._cp = self._cp / self._media['permittivity']
        #self.epsilon = self.epsilon / self.media['permittivity']   

        #Get layer coord
    def layerIndex(self, mesh):
        '''
            Gets the first and last index of the layer
        '''
        dt=mesh.steps()

        #get index of start and end position of layer in the mesh. Check a small interval of [pos-dt/2,pos+dt/2]
        indexStart = np.where((mesh.pos>=(self.pos-dt/2.0)) & (mesh.pos<(self.pos+dt/2.0)) )        
        indexEnd = np.where((mesh.pos>=((self.pos+self.width)-dt/2.0)) & (mesh.pos<((self.pos+self.width)+dt/2.0)) )        
 
        print(indexStart)
        return indexStart,indexEnd
    
    def _calcDispersionVar(self,dt,ap,cp):
        '''
            Computes kp and bp needed to calculate E and Jp in dispersive media
        '''
        dem=(1-ap*dt/2.0)
        self._kp = (1+ap*dt/2.0)/dem
        self._bp =  sp.epsilon_0*cp*dt/dem
        return self._kp,self._bp

