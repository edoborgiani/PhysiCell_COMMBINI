from cairosvg import svg2png
import glob
import os

os.chdir("/home/researcher/Physicell_COMMBINI/output")

files=glob.glob('./SVG files/Cellular level/ImmuneF/immunefluorescence*.svg')
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
    out = "./PNG files/Cellular level/ImmuneF/immunefluorescence" + sx + ".png"
    svg2png(url=inp, write_to=out)

x=0

files=glob.glob('./SVG files/Cellular level/ImmuneF/zoom_immunefluorescence*.svg')
files.sort()

for f in files:
    sx=str(x)
    if x<100:  
        sx="0"+sx
    if x<10:
        sx="0"+sx
    x+=1
    inp = f
    out = "./PNG files/Cellular level/ImmuneF/zoom_immunefluorescence" + sx + ".png"
    svg2png(url=inp, write_to=out)