import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv
from globalz import global_variables
plt.rcParams.update({'font.size': 20})

glbl = global_variables()

f_0 = 15e9
B_0 = 0.5
w_0 = f_0*2*np.pi

data_set = open("M24_w_gradients.csv",'r')
file = csv.DictReader(data_set) 

x_pos = [] 
y_pos = []
ele = []
nu = []

for row in file:
    x_pos.append(float(row['x']))
    y_pos.append(float(row['y']))
    ele.append(float(row['nD_e-']))
    nu.append(0)#float(row['nu'])) 

X_ratio = []
Y_ratio = []
for i in range(len(x_pos)):
    X_val = glbl.plasma_frequency(ele[i])/w_0
    X_ratio.append(X_val)
    Y_ratio.append(nu[i]/w_0)
print(glbl.gyrofrequency(B_0)/w_0)

plt.scatter(x_pos, y_pos, c=X_ratio,s=10, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e-4,vmax=1e2))
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.colorbar(label='X Ratio',format='%.0e') 
plt.show()
'''

plt.scatter(x_pos, y_pos, c=Y_ratio,s=10, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e-3,vmax=1e1))
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.colorbar(label='Y Ratio',format='%.0e') 
plt.show()


resonance = np.zeros_like(x_pos)
for i in range(len(x_pos)):
    X = glbl.plasma_frequency(ele[i])/w_0
    Y = nu[i]/w_0
    Z = glbl.gyrofrequency(B_0)
    denom = ((1-X)**2+Y**2-Z*(1-X**2))**2+(Y*(1-X**2)**2+Y**2+Z)**2
    resonance[i] = np.abs(denom)
plt.scatter(x_pos, y_pos, c=resonance,s=10, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e5,vmax=1e10))
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.colorbar(label='Denominator',format='%.0e') 
plt.show()

esp_m1 = np.zeros_like(x_pos)
for i in range(len(x_pos)):
    X = glbl.plasma_frequency(ele[i])/w_0
    Y = nu[i]/w_0
    Z = glbl.gyrofrequency(B_0)
    esp = X**2/(1-1j*Y-Z**2/(1-1j*Y-X**2))
    esp_m1[i] = np.abs(esp)
plt.scatter(x_pos, y_pos, c=esp_m1,s=10, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e-15,vmax=1e0))
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.colorbar(label='cutoff zero',format='%.0e') 
plt.show()
'''