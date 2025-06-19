
import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
from Input_Data.Ray_tracing_data_management import import_data_refractive_index
from Ray_Tracing_Alg.Ray_Tracing_EXE import execute_RTA 
from ISSM_Alg.ISSM_EXE import execute_ISSM
from data_plotting.plotting_program_results import result_plotting
#from Friis_components import Antenna_Parameters
#from globalz import global_variables as glbl
from INPUTS import input_data as i_d
plt.rcParams.update({'font.size': 12})
print('Running program...')

idri = import_data_refractive_index(i_d.cfd_data_file)
idri.input_variable_parameters(i_d.Frequency,i_d.B_field,m=0.176327,b=0.0213131357)

#idri.plot_cfd_data()
#idri.n_gradient()

ray_data,n_real,n_imag = execute_RTA(i_d.X_initial,i_d.Y_initial,i_d.Artificial_Antenna,idri,
angle_range=i_d.angle_range,rays=i_d.rays,refinement=i_d.refinement,step=i_d.raytracing_step_size)
#print(ray_data)
for i in range(i_d.rays):
    ray_r = ray_data[i]
    x_dat = np.array([ray[0] for ray in ray_r])
    y_dat = np.array([ray[1] for ray in ray_r])
    plt.plot(x_dat,y_dat,'r')
idri.ray_plot_background()
plt.show()

#ray_data = zip(x_dat,y_dat)

print('Calculating scattering matrix method...')
Trans_set = []
Refl_set = []

#data_set = pd.read_csv('Ray_Data.csv',low_memory=False)
for i in range(i_d.rays):
    ray_n = ray_data[i]
    #print(len(ray_n))
    n_real_set = n_real[i]
    n_imag_set = n_imag[i]
    A,D,nil,nil1 = execute_ISSM(ray_n,n_real_set,n_imag_set,i_d.Artificial_Antenna,i_d.Frequency)
    Refl_set.append(A)
    Trans_set.append(D)
print(Refl_set)
print(Trans_set)

print('Scattering matrix method solved')

print('###########################################')
print('### Transmission')
trs = np.average(Trans_set)*100
print('### '+str(trs)+'%')
print('###########################################')
print()
print('###########################################')
print('### Reflection')
ref = np.average(Refl_set)*100
print('### '+str(ref)+'%')
print('###########################################')
print()
print('###########################################')
print('### Absorption')
abso = (1-np.average(Refl_set)-np.average(Trans_set))*100
print('### '+str(abso)+'%')
print('###########################################')

print('Plotting results...')
result_plotting(ray_data,idri,trs,ref,abso)

print('Program complete')