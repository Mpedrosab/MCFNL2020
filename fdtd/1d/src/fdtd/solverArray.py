
import math
import numpy as np
import scipy.constants as sp
import copy
import time
#import matplotlib.pyplot as plt

L = 0 # Lower
U = 1 # Upper

class Fields: 
    def __init__(self, e, h):
        self.e = e
        self.h = h

    def get(self):
        return (self.e, self.h)

class ComplexField: 
    '''
        For dispersive media
    '''
    def __init__(self, Jp_old):
        self.Jp_old= Jp_old
    def get(self):
        return (self.Jp_old)        #Define for the whole mesh but only use it
                                    # where the dispersive layers is placed


class Solver:
    
    _timeStepPrint = 100

    def __init__(self, mesh, options, probes,sources, initialCond=[{"type": "none"}],dispLayer=None):
        self.options = options
        
        self._mesh = copy.deepcopy(mesh)
        self._initialCond = copy.deepcopy(initialCond)
        self._dispLayer = None          #need to check along the code if there is a layer
        #Get dispersivee media properties in case it is defined
        if dispLayer is not None:

            self._dispLayer=copy.deepcopy(dispLayer)
            self._epsilon=self._dispLayer.epsilon
            self._layerIndices = self._dispLayer.indices


            #Copy outside since it is needed in every loop. Speed up the code
            self._ap=np.zeros(( mesh.pos.size, len(self._dispLayer.ap)))
            self._cp=np.zeros(( mesh.pos.size, len(self._dispLayer.ap)))
            self._ap[self._dispLayer.indices,:]= self._dispLayer.ap
            self._cp[self._dispLayer.indices,:]= self._dispLayer.cp 
             #Changed?
            self.oldJp = ComplexField(Jp_old = np.zeros(( mesh.pos.size, len(self._dispLayer.ap))) )      #Takes the size of the layer, not the full grid
          #  self.oldJp = ComplexField(Jp_old = np.zeros(( mesh.pos.size, len(self._dispLayer.ap))) )      #Takes the size of the layer, not the full grid

        self._probes = copy.deepcopy(probes)
        for p in self._probes:
            box = self._mesh.elemIdToBox(p["elemId"])
            box = self._mesh.snap(box)
            ids = self._mesh.toIds(box)
            Nx = abs(ids)

            p["mesh"] = {"origin": box[L], "steps": abs(box[U]-box[L]) / Nx}
            p["indices"] = ids
            p["time"]   = [0.0]
            
            for initial in self._initialCond:
                if initial["type"]=="none":
                    values=np.zeros( mesh.pos.size )
                    p["values"] = [np.zeros((1,Nx[1]))]          
                elif ( initial["type"] == "gaussian"):
                    position=self._mesh.pos
                    #print(source["index"])  #Lugar del pico
                    values=Solver.movingGaussian(position, 0, \
                       sp.speed_of_light,initial["peakPosition"],\
                       initial["gaussianAmplitude"], \
                       initial["gaussianSpread"] )  
                    p["values"]= [values[ids[0]:ids[1]]]
                        #plt.plot(position,eNew)
                else:
                    raise ValueError(\
                    "Invalid initial condition type: " + initial["type"] )


        self._sources = copy.deepcopy(sources)
        for source in self._sources:
            box = self._mesh.elemIdToBox(source["elemId"])
            ids = mesh.toIds(box)
            source["index"] = ids

        self.old = Fields(e = values.copy(),
                          h = np.zeros( mesh.pos.size-1 ) )


    def solve(self, finalTime):
        tic = time.time()
        t = 0.0
        dt = self._dt()
        print ('time step:',dt)
        numberOfTimeSteps = int(finalTime / dt)

        if self._dispLayer is not None:
            self._kp,self._bp=self._dispLayer ._calcDispersionVar(dt,self._ap,self._cp)

        for n in range(numberOfTimeSteps):
            self._updateE(t, dt)
            t += dt/2.0
            self._updateH(t, dt)
            t += dt/2.0
            self._updateProbes(t)
    
            if n % self._timeStepPrint == 0 or n+1 == numberOfTimeSteps:
                remaining = (time.time() - tic) * \
                    (numberOfTimeSteps-n) / (n+1)
                min = math.floor(remaining / 60.0)
                sec = remaining % 60.0
                print("    Step: %6d of %6d. Remaining: %2.0f:%02.0f"% (n, \
                    numberOfTimeSteps-1, min, sec))
        
        print("    CPU Time: %f [s]" % (time.time() - tic))


    def _dt(self):
        return self.options["cfl"] * self._mesh.steps() / sp.speed_of_light  

    def timeStep(self):
        return self._dt()

    def getProbes(self):
        res = self._probes
        return res

    def _updateE(self, t, dt):
        (e, h) = self.old.get()
        eNew = np.zeros( self.old.e.shape )

        if self._dispLayer is None:
            cE = dt / sp.epsilon_0 / self._mesh.steps()
            eNew[1:-1] = e[1:-1] +    cE * (h[1:] - h[:-1])

        else:
            Jp_old = self.oldJp.get()
            JpNew = np.zeros( Jp_old.shape )
            cE = (2*dt/((2*self._epsilon*sp.epsilon_0+np.sum(2*np.real(self._bp),1))*self._mesh.steps()))
            #Term multiplying e[1:-1] is 1 since conductivity=0
            #Need to add an extra term to indices for dh/dx
            #indH= np.concatenate(([self._layerIndices[0]-1],self._layerIndices,))
            indH = self._layerIndices[:-1]
            #Changed?
            eNew[1:-1] = e[1:-1] + cE[1:-1] * ((h[1:] \
               - h[:-1])-np.real(np.sum((1+self._kp[1:-1,:])*Jp_old[1:-1,:],1)))
           # eNew[1:-1] = e[1:-1] + cE2 * ((h[1:] - h[:-1])-np.real(np.sum((1+self._kp)*Jp_old[1:-1],1)))


        #add mur conditions to layer
        # Boundary conditions
        for lu in range(2):
            if lu == 0:
                pos = 0
            else:
                pos = -1
            if self._mesh.bounds[lu] == "pec":
                eNew[pos] = 0.0
            elif self._mesh.bounds[lu] == 'pmc':
                eNew[pos] = e[pos] + 2*cE*(h[pos] if pos == 0 else -h[pos])
            elif self._mesh.bounds[lu] == 'mur':
                if pos == 0:
                    eNew[0] =  e[ 1]+(sp.speed_of_light*dt-self._mesh.steps())* (eNew[ 1]-e[ 0]) / (sp.speed_of_light*dt+self._mesh.steps())
                else:
                    eNew[-1] = e[-2]+(sp.speed_of_light*dt-self._mesh.steps())* (eNew[-2]-e[ 1]) / (sp.speed_of_light*dt+self._mesh.steps())
            else:
                raise ValueError("Unrecognized boundary type")

        # Source terms
        for source in self._sources:
            if source["type"] == "dipole":
                magnitude = source["magnitude"]
                if magnitude["type"] == "gaussian":
                    eNew[source["index"]] += Solver._gaussian(t, \
                        magnitude["gaussianDelay"], \
                        magnitude["gaussianSpread"] )       
                else:
                    raise ValueError(\
                    "Invalid source magnitude type: " + magnitude["type"])

            elif source["type"] == "none":
                continue
            else:
                raise ValueError("Invalid source type: " + source["type"])
        
        if self._dispLayer is not None:

            #Get new JpNew
            for i in range(0,np.shape(Jp_old)[1]):
                
                #Changed?
               JpNew[1:-1,i] = self._kp[1:-1,i]*Jp_old[1:-1,i]+ self._bp[1:-1,i] * (eNew[1:-1] - e[1:-1]) / dt
               # JpNew[1:-1,i] = self._kp[i]*Jp_old[1:-1,i]+ self._bp[i] * (eNew[1:-1] - e[1:-1]) / dt
            Jp_old[:] = JpNew[:]
        e[:] = eNew[:]
        
        
    def _updateH(self, t, dt):      
        hNew = np.zeros( self.old.h.shape )
        (e, h) = self.old.get()
        cH = dt / sp.mu_0 / self._mesh.steps()
        hNew[:] = h[:] + cH * (e[1:] - e[:-1])
        h[:] = hNew[:]
            
    def _updateProbes(self, t):
        for p in self._probes:
            if "samplingPeriod" not in p or \
               "samplingPeriod" in p and \
               (t/p["samplingPeriod"] >= len(p["time"])):
                p["time"].append(t)
                ids = p["indices"]
                values = np.zeros(ids[U]-ids[L])
                values[:] = self.old.e[ ids[0]:ids[1] ]
                p["values"].append(values)

    def _calcDispersionVar(self,dt,ap,cp):

        dem=(1-ap*dt/2.0)
        self._kp = (1+ap*dt/2.0)/dem
        self._bp =  sp.epsilon_0*cp*dt/dem
        return self._kp,self._bp

    @staticmethod
    def _gaussian(x, delay, spread):
        return np.exp( - ((x-delay)**2 / (2*spread**2)) )
    
    def movingGaussian(x,t,c,center,A,spread):
        return A*np.exp(-(((x-center)-c*t)**2 /(2*spread**2)))
