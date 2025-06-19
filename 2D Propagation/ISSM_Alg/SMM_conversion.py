import numpy as np
from globalz import global_variables
glbl = global_variables()

def convert_1D_SMM(ray_list,antenna):
    x_aa_1,y_aa_1 = antenna[0]
    x_aa_2,y_aa_2 = antenna[1]
    x_min = min([x_aa_1,x_aa_2])
    x_max = max([x_aa_1,x_aa_2])
    ray_primary = ray_list
    # print(ray_primary)
    Ray_R = [0]
    length = 0
    for i in range(0,len(ray_primary)-1):
        #print(ray_primary[i])
        P1_x1,P1_y1 = ray_primary[i]
        P1_xp1,P1_yp1 = ray_primary[i+1]
        length += np.sqrt((P1_x1-P1_xp1)**2+(P1_y1-P1_yp1)**2)
        Ray_R.append(length)
    x_fin, y_fin =ray_primary[-1]
    #print(len(Ray_R))
    #print(Ray_R[-1])
    if x_min<x_fin<x_max and y_fin>0.0:
        
        ret = 0
        #print(ret)
    else:
        ret = 1
    return Ray_R, ret

def convert_k_SMM(n_real,n_imag,freq):
    k_list = []
    for i in range(len(n_real)):
        n = n_real[i]-1j*n_imag[i]
        k_list.append(glbl.propagation_constant(freq,n))
    return k_list