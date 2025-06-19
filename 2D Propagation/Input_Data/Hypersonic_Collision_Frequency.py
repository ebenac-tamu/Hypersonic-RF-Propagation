from globalz import global_variables
import numpy as np
glbl = global_variables()

def integral(t_data,cs_data,t_e,N_i):
    result = 0
    for i in range(len(t_data)-1):
        v_th = np.sqrt(2*t_data*1.6e-19/glbl.m_e)
        #print(v_th)
        col_val = cs_data[i]*N_i*v_th
        #print(col_val)
        col_val_p1 = cs_data[i+1]*N_i*v_th
        esp_i = t_data[i]/(t_e/11600)
        #print(esp_i)
        esp_ip1 = t_data[i+1]/(t_e/11600)
        d_esp = esp_ip1-esp_i
        func_i = col_val*esp_i**(3/2)*np.exp(-1*esp_i)
        func_ip1 = col_val_p1*esp_ip1**(3/2)*np.exp(-1*esp_ip1)
        result += d_esp*(func_i+func_ip1)/2
    return result

def calculate_collisions(Temp_vib,number_density):
    # number_density is a list of number densities in the order [N,N2,O,O2,NO]
    Temp_vib*=11600
    neutral_data = ['N','N2','O','O2','NO']
    data_name = '_elasticCollisions.txt'
    folder_name = 'Input_data/'
    collision_frequency = 0
    for i in range(5):
        num_dens = number_density[i]
        Te = Temp_vib
        temp = []
        cs = []
        file_name = folder_name+neutral_data[i]+data_name
        with open(file_name, 'r') as file:
            for line in file:
                # Split each line into columns based on whitespace
                parts = line.split()
                
                # Convert the parts to float and append to respective lists
                temp.append(float(parts[0]))
                cs.append(float(parts[1]))
        #print(temp)
        #print(cs)
        product1 = 4/(3*np.sqrt(np.pi))
        product2 = integral(temp,cs,Te,num_dens)
        #print(product2)
        collision_frequency += product1*product2
    #print(collision_frequency)
    return collision_frequency