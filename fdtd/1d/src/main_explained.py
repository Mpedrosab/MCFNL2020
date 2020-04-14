# -*- coding: utf-8 -*-

#Le mete info desde archivo de entrada en formato json

import json
import argparse
import os.path
import sys

from fdtd.mesh import Mesh
from fdtd.solver import Solver
from fdtd.viewer import Animator

print("=== Python FDTD 1D")

parser = argparse.ArgumentParser(description='Python FDTD 1D')
parser.add_argument('-i', '--input', nargs=1, type=str)
args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()

inputFilename = ''.join(args.input).strip()
print("--- Reading file: %s"%(inputFilename))
data = json.load(open(inputFilename))               #Lee los datos en formato json

print('--- Initializing mesh')
#Utiliza la clase Mesh para inicializar la malla (ver la clase para info)
mesh = Mesh(data["coordinates"], data["elements"], data["grid"])    #Initialize grid
#Pega: me tengo que aprender qué significa cada etiqueta

print('--- Initializing solver')
solver = Solver(mesh, data["options"], data["probes"], data["sources"])         #Todo el meollo, c+ómo se actualiza los campos, cambios...
                                                                                                                                                                                                                        

print('--- Solving')
solver.solve(data["options"]["finalTime"])

print('--- Visualizing')
Animator(mesh, solver.getProbes()[0], skipFrames=200)

print('=== Program finished')