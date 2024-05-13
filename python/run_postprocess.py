# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:18:04 2022

@author: u0137011
"""

import os
import sys
import numpy as np

py_folder="$HOME/host_home/Physicell_COMMBINI/python"
pv_folder="$ROOT/ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64/bin"

list_cmd=["python VTK_generator.py","python SVG_to_PNG.py","python SVG_to_PNG_IF.py","python GIF_maker.py",pv_folder + "/pvpython VTKtoSVG_pvpy.py","python SVGmerger.py","python SVG_to_PNG_microenv.py","python GIF_maker_microenv.py","python cell_counter.py"]

i=0
sys.stdout.write("[%-20s] %d%%" % ('='*int(i*20/len(list_cmd)), i*100/len(list_cmd)))
sys.stdout.flush()
i+=1

for c in list_cmd:
    os.system(c)
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d%%" % ('='*int(i*20/len(list_cmd)), i*100/len(list_cmd)))
    i+=1
    sys.stdout.flush()

sys.stdout.write('\n')
sys.stdout.flush()