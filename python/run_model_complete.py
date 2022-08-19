# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:18:04 2022

@author: u0137011
"""

import os

py_folder="$HOME/host_home/Physicell_COMMBINI/python"
pv_folder="$ROOT/ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64/bin"

os.system("python run_model.py")

os.system("python VTK_generator.py")

os.system("python SVG_to_PNG.py")

os.system("python SVG_to_PNG_IF.py")

os.system("python GIF_maker.py")

os.system(pv_folder + "/pvpython VTKtoSVG_pvpy.py")