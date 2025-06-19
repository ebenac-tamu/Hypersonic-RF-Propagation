import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv
from globalz import global_variables
plt.rcParams.update({'font.size': 20})

glbl = global_variables()

B_n = lambda X,Y,Z: (X**4-2*X**2+Z**2*X**2-Z**2+Y**2+1)-1j*(X**4*Y-2*X**2*Y+Y**3+Y+Z**2*Y)
A_n = lambda X,Y,Z: X**2*(X**4-2*X**2+Y**2+1)

dBdNe = lambda X,Y,Z: (4*X**3-4*X+Z**2*X)-1j*(4*X**3*Y-4*X*Y)
dBdNu = lambda X,Y,Z: 2*Y-1j*(X**4-2*X**2+3*Y**2+Z**2+1)
dAdNe = lambda X,Y,Z: 6*X**5-8*X**3+2*(1+Y**2)*X
dAdNu = lambda X,Y,Z: 2*X**2*Y

n2_AH = lambda X,Y,Z: 1-X**2/(1-1j*Y-Z**2/(1-1j*Y-X**2))

def grad_n(Ne,Nu,f0,B_mag,grad_Wp,grad_Nu,grad_n0):
    Wp = glbl.plasma_frequency(Ne)
    Wc = glbl.gyrofrequency(B_mag)
    omega = 2*np.pi*f0
    X0 = Wp/(omega)
    Y0 = Nu/(omega)
    Z0 = Wc/(omega)
    Wp_x = grad_Wp[0]
    Wp_y = grad_Wp[1]
    Wp_z = grad_Wp[2]
    Nu_x = grad_Nu[0]
    Nu_y = grad_Nu[1]
    Nu_z = grad_Nu[2]
    n0_x = grad_n0[0]
    n0_y = grad_n0[1]
    n0_z = grad_n0[2]
    A_0 = A_n(X0,Y0,Z0)
    B_0 = B_n(X0,Y0,Z0)
    A_Ne = dAdNe(X0,Y0,Z0)
    B_Ne = dBdNe(X0,Y0,Z0)
    A_Nu = dAdNu(X0,Y0,Z0)
    B_Nu = dBdNu(X0,Y0,Z0)
    n2 = n2_AH(X0,Y0,Z0)
    
    #print(B_0)

    dndWp = (A_Ne*B_0-A_0*B_Ne)/(2*omega*B_0**2) / np.sqrt(n2)
    dndNu = (A_Nu*B_0-A_0*B_Nu)/(2*omega*B_0**2) / np.sqrt(n2)
    
    grad_n = np.array([0,0,0])
    grad_n[0] = np.real(n0_x + np.real(dndWp*Wp_x) + np.real(dndNu*Nu_x))
    grad_n[1] = np.real(n0_y + np.real(dndWp*Wp_y) + np.real(dndNu*Nu_y))
    grad_n[2] = np.real(n0_z + dndWp*Wp_z + dndNu*Nu_z)
    #print('grad x')
    #print(np.abs(grad_n[0])**2)
    #print('grad y')
    #print(np.abs(grad_n[1])**2)
    grad_n_mag = np.sqrt(grad_n[0]**2+grad_n[1]**2)
    return grad_n,grad_n_mag,np.sqrt((np.real(dndNu*Nu_x)+n0_x)**2 + (np.real(dndNu*Nu_x)+n0_y)**2)
    

print('Extracting CFD Data...')
data_set = open("M24_w_gradients.csv",'r')
file = csv.DictReader(data_set) 
x_pos = [] 
y_pos = []

n_0 = []
grad_wp = []
del_wp_mag = []
grad_neutral = []
grad_col = [] 
ele = []
Temp_vib = []
rho_air = [] 
nu = []
for row in file:
    #print(row)
 
    x_pos.append(float(row['x']))

    y_pos.append(float(row['y']))
    grad_wp.append([float(row['Gradient_Wp_x']),float(row['Gradient_Wp_y']),float(row['Gradient_Wp_z'])])
    
    grad_col.append([float(row['Gradient_col_x']),float(row['Gradient_col_y']),float(row['Gradient_col_z'])])
    grad_neutral.append([float(row['Gradient_neutral_x']),float(row['Gradient_neutral_y']),float(row['Gradient_neutral_z'])])
    del_wp_mag.append(np.sqrt(float(row['Gradient_neutral_x'])**2 + float(row['Gradient_neutral_y'])**2))
    n_0.append(float(row['refra_neutral']))
 
    ele.append(float(row['nD_e-']))

    nu.append(float(row['nu']))   

B_0 = 0.5
f_0 = 15e9


gradient_n = []
gradient_x = []
gradient_mag = []
dndne = []
for i in range(len(x_pos)):
    #print(grad_wp[i])
    del_n,del_mag,dne = grad_n(ele[i],nu[i],f_0,B_0,grad_wp[i],grad_col[i],grad_neutral[i])
    #print(del_n)
    gradient_n.append(del_n)
    gradient_x.append(del_n[0])
    gradient_mag.append(del_mag)
    dndne.append(dne)

plt.scatter(x_pos, y_pos, c=dndne,s=10, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e-10,vmax=1e4))
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.colorbar(label='grad(n)',format='%.0e') 
plt.show()