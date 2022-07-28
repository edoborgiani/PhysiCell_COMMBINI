from cairosvg import svg2png
import os

os.chdir("C:/Workdir/Programs/PhysiCell_V.1.7.0/PhysiCell/output")

for x in range(121):
    sx=str(x)
    if x<100:  
        sx="0"+sx
    if x<10:
        sx="0"+sx
    inp = "./SVG files/Cellular level/snapshot00000" + sx + ".svg"
    out = "./PNG files/Cellular level/snapshot00000" + sx + ".png"
    svg2png(url=inp, write_to=out)

for x in range(121):
    sx=str(x)
    if x<100:  
        sx="0"+sx
    if x<10:
        sx="0"+sx
    inp = "./SVG files/Cellular level/zoom00000" + sx + ".svg"
    out = "./PNG files/Cellular level/zoom00000" + sx + ".png"
    svg2png(url=inp, write_to=out)