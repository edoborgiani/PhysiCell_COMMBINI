#%%
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 14:18:04 2022

@author: u0137011
"""

import glob
import shutil
import os

list_microenv=['Debris','IFNg','TNFa','IL10','TGFb']

for microelem in list_microenv:

    files_svg=glob.glob("../output/SVG files/Microenvironment/" + microelem + "_*.svg")
    files_svg.sort()

    #print(files_svg)

    for x in range(0,len(files_svg)):
        sx=str(x)
        if x<100:  
            sx="0"+sx
        if x<10:
            sx="0"+sx
        
        f_cell="../output/SVG files/Cellular level/snapshot00000"+sx+".svg"
        f_menv="../output/SVG files/Microenvironment/" + microelem + "_"+sx+".svg"
        f_total="../output/SVG files/" + microelem + "_" + sx + ".svg"
        f_app="./append_file2.txt"

        shutil.copyfile(f_app,f_total)

        #os.system("chmod 777 " + f_total)

        with open(f_total,"a") as finfile:
            finfile.write('\n\t<g transform="matrix(7.4685891,0,0,7.2597233,-502.10967,-311.14899)" shape-rendering="crispEdges">')
            #for line in finfile:
                #print(line)

        lincount=0
        
        with open(f_total,'a') as writefile, open(f_menv, 'r', encoding='utf-8') as readfile:
            for line in readfile:
                #print(line)
                if lincount>0:
                    writefile.write('\t'+line)
                    lincount-=1

                if line.startswith('<g>'):
                    #print(line)
                    lincount=3
                
            writefile.write('\t</g>')

        linvalid=False

        with open(f_total,'a') as writefile, open(f_cell, 'r', encoding='utf-8') as readfile:
            for line in readfile:
                #print(line)

                if line.startswith('<pa'):
                    #print(line)
                    linvalid=True

                if linvalid:
                    writefile.write('\t'+line)