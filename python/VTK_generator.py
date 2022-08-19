#%%
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:18:04 2022

@author: u0137011
"""
import os
import time
from os import path
import numpy as np
from scipy.io import loadmat  # this is the SciPy module that loads mat-files
import matplotlib.pyplot as plt
from datetime import datetime, date, time
import pandas as pd
from pyevtk.hl import gridToVTK
import glob

curpath=os.path.abspath(os.getcwd()).replace("\\","/")
print(curpath)

completed=0

os.chdir("..")

files=glob.glob('./output/new_*_microenvironment0.mat')
files.sort()
b=0

for f in files:

    sx=str(b)
    if b<100:  
        sx="0"+sx
    if b<10:
        sx="0"+sx
    b+=1

    mat = loadmat(f)  # load mat-file
    A=mat['multiscale_microenvironment']
    #print(A)
    xmin = A[0,0]
    ymin = A[1,0] 
    zmin = A[2,0] 

    #print(xmin)

    n = len(A[0])-1
    xmax = A[0,n] 
    ymax = A[1,n] 
    zmax = A[2,n] 

    #print(xmax)

    xnodes = 0
    while( A[0,xnodes] < xmax-1e-10): 
        xnodes+=1

    xnodes+=1

    X = A[0,0:xnodes+1]

    #print(xnodes)

    ynodes = 0
    while( A[1,ynodes*xnodes] < ymax-1e-10): 
        ynodes+=1

    ynodes+=1

    Y = A[1,0:xnodes*ynodes:xnodes]

    #print(ynodes)

    znodes = 0
    while( A[2,znodes*xnodes*ynodes] < zmax-1e10): 
        znodes+=1

    znodes+=1    

    Z = A[2,0:znodes*xnodes*ynodes:xnodes*ynodes]

    #print(znodes)

    temp = np.zeros((xnodes,ynodes,znodes,5))

    for n in range(4,len(A)):
        m=0
        for k in range(0,znodes):
            for j in range(0,ynodes):
                for i in range(0,xnodes):
                    temp[i,j,k,n-4]=A[n,m]
                    m+=1

    #plt.imshow(temp[...,2])
    #plt.colorbar()
    #plt.savefig("test.png")

    output_debris = np.dstack([temp[...,0]])
    output_TNF = np.dstack([temp[...,1]])
    output_TGF = np.dstack([temp[...,2]])
    output_IL10 = np.dstack([temp[...,3]])
    output_IFN = np.dstack([temp[...,4]])

    x = np.arange(0, xnodes+1)
    y = np.arange(0, ynodes+1)
    z = np.arange(0, znodes)

    gridToVTK("./output/VTK files/microenvironment_" + sx ,  x, y, z, cellData={"Debris":output_debris, "TNFa":output_TNF, "IL10":output_IL10, "TGFb":output_TGF, "IFNg":output_IFN})