# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:18:04 2022

@author: u0137011
"""

import os
import sys
import numpy as np

py_folder="$HOME/host_home/Physicell_COMMBINI/python"
output_folder="$HOME/host_home/Physicell_COMMBINI/output"
output_stable_folder="$HOME/host_home/outputs"
pv_folder="$ROOT/ParaView-5.10.1-osmesa-MPI-Linux-Python3.9-x86_64/bin"

for k in ['a','b','c','d','e']:
    os.system("python run_model_complete.py")

    os.system("cp " + output_folder + "/result.csv " + output_stable_folder + "/result_semirigid_mD_" + k + ".csv")
    for m in range(3):
        os.system("cp " + output_folder + "/'CSV files'/cells_total_0" + str(m+1) + ".csv " + output_stable_folder + "/cells_total_semirigid_mD_" + k + "0" + str(m+1) + ".csv")

