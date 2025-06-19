import matplotlib.pyplot as plt
import numpy as np

def ray_plot(ray_data,idri):
    x_plot = []
    y_plot = []
    #n_real = []
    for i in range(len(ray_data)):
        x_value,y_value = ray_data[i]
        #n_val, nul1, nul2 = idri.refractive_index_set(x_value,y_value,[0,0,0])
        x_plot.append(x_value)
        y_plot.append(y_value)
        #n_real.append(np.real(n_val))
    plt.plot(x_plot,y_plot)
    #plt.scatter(x_plot,y_plot,c=n_real,cmap='jet',s=0.5)