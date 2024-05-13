from cairosvg import svg2png
import glob
import os

os.chdir("/home/researcher/Physicell_COMMBINI/output")
#print(os.listdir())

files=glob.glob('./SVG files/Cellular level/snapshot*.svg')
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
    out = "./PNG files/Cellular level/snapshot00000" + sx + ".png"
    svg2png(url=inp, write_to=out)

x=0

files=glob.glob('./SVG files/Cellular level/zoom*.svg')
files.sort()

for f in files:
    sx=str(x)
    if x<100:  
        sx="0"+sx
    if x<10:
        sx="0"+sx
    x+=1
    inp = f
    out = "./PNG files/Cellular level/zoom00000" + sx + ".png"
    svg2png(url=inp, write_to=out)