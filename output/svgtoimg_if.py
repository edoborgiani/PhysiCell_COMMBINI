from cairosvg import svg2png
import os

os.chdir("C:/Workdir/Programs/PhysiCell_V.1.7.0/PhysiCell/output")

for x in range(20):
    sx=str(x)
    if x<10:
        sx="0"+sx
    inp = "./SVG files/Cellular level/ImmuneF/immunefluorescence" + sx + ".svg"
    out = "./PNG files/Cellular level/ImmuneF/immunefluorescence" + sx + ".png"
    svg2png(url=inp, write_to=out)

for x in range(20):
    sx=str(x)
    if x<10:
        sx="0"+sx
    inp = "./SVG files/Cellular level/ImmuneF/zoom_immunefluorescence" + sx + ".svg"
    out = "./PNG files/Cellular level/ImmuneF/zoom_immunefluorescence" + sx + ".png"
    svg2png(url=inp, write_to=out)