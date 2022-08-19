import imageio
import glob
import os

os.chdir("/root/host_home/Physicell_COMMBINI/output")

files=glob.glob('./PNG files/Cellular level/snapshot*.png')
files.sort()

final_gif=[]

for f in files:
    final_gif.append(imageio.v3.imread(f))

imageio.mimsave('./PNG files/Cellular level/animation.gif', final_gif, duration=0.1)

files=glob.glob('./PNG files/Cellular level/zoom*.png')
files.sort()

final_gif=[]

for f in files:
    final_gif.append(imageio.v3.imread(f))

imageio.mimsave('./PNG files/Cellular level/animation_zoom.gif', final_gif, duration=0.1)