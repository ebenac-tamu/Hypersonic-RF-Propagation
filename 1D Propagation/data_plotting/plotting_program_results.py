import matplotlib.pyplot as plt
from data_plotting.ray_plotting import ray_plot
from data_plotting.color_plotting import color_plot
from data_plotting.bar_graph import SMM_bar_graph

def result_plotting(rays,data_import,Transmission,Reflection,Absorption):
    x_dat,y_dat,ele = data_import.extract_data_for_plotting()
    color_plot(x_dat,y_dat,ele,log_scale=True)
    for i in range(len(rays)):
        ray_plot(rays[i],data_import)
    #plt.colorbar(label='n-1')
    plt.show()
    SMM_bar_graph(Transmission,Reflection,Absorption)
    plt.show()