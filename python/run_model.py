# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:18:04 2022

@author: u0137011
"""

import os
import time
from os import path

curpath=os.path.abspath(os.getcwd()).replace("\\","/")
print(curpath)

completed=0

os.chdir("C:/Workdir/Programs/PhysiCell_V.1.7.0/PhysiCell_2")

for n in range(1,2,4):

    os.system("make reset")
    os.system("make clean")
    os.system("make data-cleanup")
    os.system("make macrophages")
    os.system("make")
    while (completed==0):
        os.system("project_migration.exe")
        os.system("timeout /t 5")
        if path.exists("./output/CSV files/cells_final.csv"):
            completed=1
        print(completed)
    
    completed=0
        
    #os.system(resultcmd)
    os.system("timeout /t 5")


        
    