import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
from Input_Data.Ray_tracing_data_management import import_data_refractive_index
from ISSM_Alg.ISSM_EXE import execute_ISSM
#from Friis_components import Antenna_Parameters
#from globalz import global_variables as glbl
from SMM_input import input_data as i_d
import csv
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.family'] = 'Times New Roman'
print('Running program...')
#y_value = lambda x: (i_d.Y_range[1]-i_d.Y_range[0])/(i_d.X_range[1]-i_d.X_range[0])*(x-i_d.X_range[0])+i_d.Y_range[0]-0.0075
#print(y_value(i_d.X_input))
#print((i_d.Y_range[1]-i_d.Y_range[0])/(i_d.X_range[1]-i_d.X_range[0])*(-1*i_d.X_range[0])+i_d.Y_range[0])
idri = import_data_refractive_index(i_d.cfd_data_file)
idri.input_variable_parameters(i_d.Frequency,i_d.B_field,height=i_d.height)
idri.input_variable_parameters(i_d.Frequency,i_d.B_field,i_d.height)
# Uncomment here to visualize field data
###########################
#idri.plot_cfd_data()
###########################

### heights 0.012,0.014,0.016

### Point 1 is at -0.053142,0.00612861   \\\\\\\\\\\\\\\\\\\\ 0.0997,0.038904
### Point 2 is at -0.0101309998,0.0114100015
### Point 3 is at -0.02502299, 0.009581216   \\\\\\\\\\\\\\\\ 0.5008,0.10962

x_range = np.linspace(i_d.X_range[0],i_d.X_range[-1],200)
x_plot = np.linspace(i_d.X_range[0],i_d.X_range[-1],200)
y_value = lambda x: (i_d.Y_range[1]-i_d.Y_range[0])/(i_d.X_range[1]-i_d.X_range[0])*(x-i_d.X_range[0])+i_d.Y_range[0]
Y_input = y_value(i_d.X_input)

T = []
R = []
A = []

x_pos = []
freq = []

height = lambda x: (0.12-0.075)/(i_d.X_range[1]-i_d.X_range[0])*(x-i_d.X_range[0])+0.075

if i_d.run_type=='PosComp':
    plt.figure(figsize=(3.2,2.8))
    for i in range(len(x_range)):
        idri.input_variable_parameters(i_d.Frequency,i_d.B_field,height(x_range[i])) 
        n_real,n_imag,length = idri.extract_SMM_inputs(x_range[i],y_value(x_range[i]))
        #plt.plot(length,n_real,label='n_real')
        #plt.plot(length,n_imag,label='n_imag')
        #plt.legend()
        #plt.show()
        A_mag,D_mag,A_phase,D_phase,b_val,c_val,B_phase,C_phase = execute_ISSM(length,n_real,n_imag,i_d.Frequency)
        #print(D_mag)
        T.append(D_mag *100)
        R.append(A_mag *100)
        A.append(100-T[i]-R[i])
    #valid_indices = [i for i in range(len(T)) if not math.isnan(T[i])]
    #x_valid = [x_range[i] for i in valid_indices]
    #T_valid = [T[i] for i in valid_indices]
    #R_valid = [R[i] for i in valid_indices]
    #A_valid = [A[i] for i in valid_indices]
    plt.plot(x_plot,T,label='Transmission')
    plt.plot(x_plot,R,label='Reflection')
    plt.plot(x_plot,A,label='Absorption')
    plt.ylabel('Power (%)')
    plt.xlim(x_plot[0],x_plot[-1])
    plt.xlabel('X-Position (m)')
    plt.ylim(0,100)
    plt.legend(loc='upper right')
    plt.subplots_adjust(left=0.214, bottom=0.179, right=0.907, top=0.924, wspace=0.2, hspace=0.2)
    plt.show()
    header = ['X','T','R','A']
    rows = zip(i_d.X_range,T,R,A)


'''
#collision comparison
fig, axs = plt.subplots(3, sharex=True)
for j in range(0,3):
    T = []
    R = []
    A = []
    for i in range(len(i_d.B_range)):
        
        idri.input_variable_parameters(i_d.Frequency,[i_d.B_range[i],0.0,0.0],i_d.height)
            
        n_real,n_imag,length = idri.extract_SMM_inputs(i_d.X_input,y_value(i_d.X_input),ij=j)
        #plt.plot(length,n_real,label='n_real')
        #plt.plot(length,n_imag,label='n_imag')
        #plt.legend()
        #plt.show()
        A_mag,D_mag,nil1,nil2 = execute_ISSM(length,n_real,n_imag,i_d.Frequency)
        #print(D_mag)
        T.append(D_mag *100)
        R.append(A_mag *100)
        A.append(100-T[i]-R[i])
    #valid_indices = [i for i in range(len(T)) if not math.isnan(T[i])]
    #x_valid = [x_range[i] for i in valid_indices]
    #T_valid = [T[i] for i in valid_indices]
    #R_valid = [R[i] for i in valid_indices]
    #A_valid = [A[i] for i in valid_indices]
    if j==0:
        axs[0].plot(i_d.B_range,T,'b',label='Transmission')
        #plt.plot(i_d.B_range,T,'b',label='Transmission')
        axs[1].plot(i_d.B_range,R,'r',label='Reflection')
        #plt.plot(i_d.B_range,R,'r',label='Reflection')
        axs[2].plot(i_d.B_range,A,'g',label='Attenuation')
        #plt.plot(i_d.B_range,A,'g',label='Attenuation')
    if j ==1:
        axs[0].plot(i_d.B_range,T,'b--',label='Transmission, No Collisions')
        #plt.plot(i_d.B_range,T,'b--',label='Transmission, No Collisions')
        axs[1].plot(i_d.B_range,R,'r--',label='Reflection, No Collisions')
        #plt.plot(i_d.B_range,R,'r--',label='Reflection, No Collisions')
        axs[2].plot(i_d.B_range,A,'g--',label='Attenuation, Increased Altitude')
    if j ==2:
        axs[0].plot(i_d.B_range,T,'b*--',label='Transmission, Decreased Altitude')
        axs[1].plot(i_d.B_range,R,'r*--',label='Reflection, Decreased Altitude')
        axs[2].plot(i_d.B_range,A,'g*--',label='Attenuation, Decreased Altitude')
# Set labels and limits

for ax in axs:
    ax.set_ylabel('Power (%)')
    ax.set_xlim(i_d.B_range[0], i_d.B_range[-1])
    ax.set_ylim(0, 100)
    ax.legend()
axs[2].set_xlabel('Magnetic Field (T)')
#plt.ylabel('Power (%)')
#plt.xlim(i_d.B_range[0], i_d.B_range[-1])
#plt.ylim(0, 100)
#plt.legend()
#plt.xlabel('Magnetic Field (T)')
plt.show()
'''
'''

x_val = [-0.04,-0.025,-0.01]
for i in x_val:
    X,Y,length = idri.Ratio_data(i,y_value(i))

    plt.plot(length,X,label='X at '+str(round(0.07+i,3))+' m')
    plt.plot(length,Y,label='Y at '+str(round(0.07+i,3))+' m')
#plt.ylabel('Power (%)')
plt.xlabel('length (m)')
#plt.xlim(0,i_d.B_range[-1])
#plt.ylim(0,100)
plt.legend()
plt.show()
'''

if i_d.run_type=='FComp':
    plt.figure(figsize=(3.2,2.8))
    for j in i_d.frequency_range:
        T = []
        R = []
        R_phase = []
        A = []
        for i in range(len(i_d.B_range)):
            
            idri.input_variable_parameters(j,[i_d.B_range[i],0.0,0.0],i_d.height)
                
            n_real,n_imag,length = idri.extract_SMM_inputs(i_d.X_input,i_d.Y_input)
            #print(n_real)
            #plt.plot(length,n_real,label='$n_{real}$')
            #plt.plot(length,n_imag,label='$n_{imag}$')
            #plt.legend()
            #plt.show()
            A_mag,D_mag,A_phase,D_phase,b_val,c_val,B_phase,C_phase = execute_ISSM(length,n_real,n_imag,j)
            #print(D_mag)
            T.append(D_mag *100)
            R.append(A_mag *100)
            A.append(100-T[i]-R[i])
            R_phase.append(A_phase)
            
        #valid_indices = [i for i in range(len(T)) if not math.isnan(T[i])]
        #x_valid = [x_range[i] for i in valid_indices]
        #T_valid = [T[i] for i in valid_indices]
        #R_valid = [R[i] for i in valid_indices]
        #A_valid = [A[i] for i in valid_indices]

        plt.plot(i_d.B_range,T,label=str(round(j/1e9,1))+'GHz')
        #plt.plot(i_d.B_range,R_phase,label=str(round(j/1e9,1))+'GHz')
        #plt.plot(i_d.B_range,R,'r',label='Reflection')
        #plt.plot(i_d.B_range,A,'g',label='Absorption')
    #plt.ylabel('Phase')
    plt.ylabel('Power (%)')
    plt.xlabel('Magnetic Field (T)')
    plt.xlim(i_d.B_range[0],i_d.B_range[-1])
    plt.ylim(0,100)
    #plt.ylim(-np.pi/2,np.pi/2)
    plt.legend(loc='lower right')
    plt.subplots_adjust(left=0.214, bottom=0.179, right=0.907, top=0.924, wspace=0.2, hspace=0.2)
    plt.show()
    header = ['F','T','R','A']
    rows = zip(i_d.frequency_range,T,R,A)


if i_d.run_type=='BComp':
    T = []
    T_phase = []
    R = []
    R_phase = []
    A = []
    for i in range(len(i_d.B_range)):
        
        idri.input_variable_parameters(i_d.Frequency,[i_d.B_range[i],0.0,0.0],i_d.height)
            
        n_real,n_imag,length = idri.extract_SMM_inputs(i_d.X_input,i_d.Y_input)
        '''
        if i==0:
            plt.plot(length,n_real,label='n_real')
            plt.plot(length,n_imag,label='n_imag')
            plt.legend()
            plt.show()
        '''
        A_mag,D_mag,A_phase,D_phase,b_val,c_val,B_phase,C_phase = execute_ISSM(length,n_real,n_imag,i_d.Frequency)
        #print(D_mag)
        T.append(D_mag *100)
        R.append(A_mag *100)
        T_phase.append(D_phase)
        R_phase.append(A_phase)
        A.append(100-T[i]-R[i])
        B_mag = np.abs(b_val)**2
        C_mag = np.abs(c_val)**2
        #C_mag_real = np.real(c_val)
        #C_mag_imag = np.imag(c_val)
        plt.plot(length,B_mag,label='B='+str(round(i_d.B_range[i],1))+' T (B)')
        #plt.plot(length,C_phase,label='B='+str(i_d.B_range[i])+' T (C)')
    #plt.ylabel('Phase')
    plt.xlabel('Position (m)')
    plt.legend()
    plt.show()
        
    #valid_indices = [i for i in range(len(T)) if not math.isnan(T[i])]
    #x_valid = [x_range[i] for i in valid_indices]
    #T_valid = [T[i] for i in valid_indices]
    #R_valid = [R[i] for i in valid_indices]
    #A_valid = [A[i] for i in valid_indices]

    plt.plot(i_d.B_range,T,label='Transmission')
    plt.plot(i_d.B_range,R,'r',label='Reflection')
    plt.plot(i_d.B_range,A,'g',label='Absorption')
    plt.ylabel('Power (%)')
    plt.xlabel('Magnetic Field (T)')
    plt.xlim(i_d.B_range[0],i_d.B_range[-1])
    plt.ylim(0,100)
    plt.legend()
    plt.show()
    plt.plot(i_d.B_range,T_phase,label='Transmission Phase')
    plt.plot(i_d.B_range,R_phase,'r',label='Reflection Phase')
    plt.ylabel('Phase')
    plt.xlabel('Magnetic Field (T)')
    plt.xlim(i_d.B_range[0],i_d.B_range[-1])
    plt.ylim(-np.pi/2,np.pi/2)
    plt.legend()
    plt.show()
    header = ['B','T','R','A']
    rows = zip(i_d.B_range,T,R,A)
    

with open("output.csv",'w',newline="") as file:
    writer = csv.writer(file)
    writer.writerow(header)
    writer.writerows(rows)
print("CSV file has been written successfully.")