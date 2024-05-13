import glob
import os
import scipy.io
import array as arr
import csv

os.chdir("/home/researcher/Physicell_COMMBINI/output")

files=glob.glob('./XML and MAT files/new_output*_cells_physicell.mat')
files.sort()

csv_tot=[]

for i in files:
    mat = scipy.io.loadmat(i)
    #mat_2 = [d for d in if d[2] > 100.0]
    #print(str(mat['cells'][2]>100.0))
    #print(mat)

    #cell_array=arr.array('f',mat['cells'][5]*(mat['cells'][4]>10.0)*(mat['cells'][1]>0.0)*(mat['cells'][2]>-500.0)*(mat['cells'][2]<500.0))
    cell_array=arr.array('f',mat['cells'][5])
    #print(cell_array)
    m0_count=cell_array.count(1)
    m1_count=cell_array.count(2)
    m2_count=cell_array.count(3)
    mn_count=cell_array.count(4)
    csv_tot.append([m0_count,m1_count,m2_count,mn_count])

fields = ['M0', 'M1', 'M2', 'PMN'] 
  
with open('result.csv', 'w') as f:
      
    # using csv.writer method from CSV package
    write = csv.writer(f)
      
    write.writerow(fields)
    write.writerows(csv_tot)


