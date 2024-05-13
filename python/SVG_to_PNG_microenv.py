from cairosvg import svg2png
import glob
import os

os.chdir("/home/researcher/Physicell_COMMBINI/output")

list_microenv=['Debris','IFNg','TNFa','IL10','TGFb']

for microenv in list_microenv:

    files=glob.glob('./SVG files/'+ microenv + '*.svg')
    files.sort()

    x=0

    for f in files:
        sx=str(x)
        if x<100:  
            sx="0"+sx
        if x<10:
            sx="0"+sx
        x+=1
        inp = f
        out = "./PNG files/Microenvironment/" + microenv + "_" + sx + ".png"
        svg2png(url=inp, write_to=out)
