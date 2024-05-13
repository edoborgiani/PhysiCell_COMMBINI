import imageio
import glob
import os

os.chdir("/home/researcher/Physicell_COMMBINI/output")

list_microenv=['Debris','IFNg','TNFa','IL10','TGFb']

for menv in list_microenv:
    files=glob.glob('./PNG files/Microenvironment/' + menv + '_*.png')
    files.sort()

    final_gif=[]

    for f in files:
        final_gif.append(imageio.v3.imread(f))

    imageio.mimsave('./PNG files/' + menv + '_animation.gif', final_gif, duration=0.1)