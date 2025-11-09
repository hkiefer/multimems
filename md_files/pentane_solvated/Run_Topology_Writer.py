import os 
import math
import shutil

gmx_file = open('topol.top','w')

# Choose the force field:
gmx_file.write(''.join(['#include "gromos54a7.ff/forcefield.itp"','\n']))

# Choose the water models:
gmx_file.write(''.join(['#include "spce_custom.itp"','\n']))

# Choose additional molecules:
gmx_file.write(''.join(['#include "pentane.itp"','\n']))

gmx_file.write(''.join(['[system]','\n']))
gmx_file.write(''.join(['glycerol in water','\n\n']))

gmx_file.write(''.join(['[ molecules ]','\n']))
gmx_file.write(''.join(['; Compound        #mols','\n']))

print('Topology file written...')




