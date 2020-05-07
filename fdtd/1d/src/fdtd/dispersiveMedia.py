import numpy as np
import copy
import math
import scipy.constants as sp

class DispersiveMedia:
    def __init__(self,mesh,wall):
        #Copy outside since it is needed in every loop. Speed up the code
        self._media=copy.deepcopy(wall)
        self.ap=np.empty(len(self._media["ap"]),dtype=complex)
        self.cp=np.empty(len(self._media["cp"]),dtype=complex)
        self.epsilon=self._media['permittivity'] * sp.epsilon_0
        for i in range(0,len(self._media["ap"])):
            self.ap[i]=complex(self._media["ap"][i])
            self.cp[i]=complex(self._media["cp"][i])
        
        self.cp=self.cp * sp.epsilon_0
    #WATCH OUT! Normalize permittivity0
        #self._ap = self._ap / self._media['permittivity']
        #self._cp = self._cp / self._media['permittivity']
        #self.epsilon = self.epsilon / self.media['permittivity']   

        #Get coord of wall

