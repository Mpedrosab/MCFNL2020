
import json
import argparse
import os.path
import sys

from fdtd.mesh import Mesh
from fdtd.solver import Solver
from fdtd.viewer import Animator
from fdtd.comparison import AnalyticComp
from fdtd.dispersiveMedia import DispersiveMedia
from measure.Transmittance import MeasureTransmittance, AnalyticTransmittance, PlotTransmittance
print("=== Python FDTD 1D")

'''
parser = argparse.ArgumentParser(description='Python FDTD 1D')
parser.add_argument('-i', '--input', nargs=1, type=str)
args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

inputFilename = ''.join(args.input).strip()
'''
inputFilename='..\\tests\\cavity_dispersive_test_NoComplex.json'
print("--- Reading file: %s"%(inputFilename))
data = json.load(open(inputFilename))



print('--- Initializing mesh')
mesh = Mesh(data["coordinates"], data["elements"], data["grid"])
#layer = None
layer = DispersiveMedia(mesh,data["dispersiveLayers"])
#indeces = layer.layerIndices(mesh)
print('--- Initializing solver')

solver = Solver(mesh, data["options"], data["probes"], data["sources"],dispLayer = layer)

print('--- Solving')
solver.solve(data["options"]["finalTime"])

print('--- Visualizing')
solNum=solver.getProbes()[0]

#Measure transmittance
transmittance = MeasureTransmittance(layer,solNum)
freq, transfft = transmittance.AmplVsFreq()
PlotTransmittance(freq, transfft)


#Analytical transmittance
transmittanceReal = AnalyticTransmittance(layer)
transReal = transmittanceReal.T(freq)
PlotTransmittance(freq, transReal)

Animator(mesh, solNum,layer=layer)
#%%
'''
print('--- Comparison with analytical solution')
comparison=AnalyticComp(mesh, solNum, data["initialCond"])
solReal=comparison.AnalyticalSol(comparison.gridE,comparison.probeTime)
err=comparison.L2Error(solReal)
comparison.AnimatorTogether(solNum,solReal)
comparison.PrintErr(solNum,err)

print('=== Program finished')
'''