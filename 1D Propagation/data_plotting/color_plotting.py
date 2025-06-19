import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def color_plot(x_points,y_points,flow_data,log_scale=False,colorbar_title=''):
    if log_scale ==True:
        min_val = min(flow_data)
        max_val = max(flow_data)
        plt.scatter(x_points, y_points, c=flow_data,s=100, cmap='gray_r', alpha=0.5, norm=LogNorm(vmin=1e12,vmax=1e22))
    else:
        plt.scatter(x_points, y_points, c=flow_data,s=100, cmap='viridis', alpha=0.5)
    plt.colorbar(label=colorbar_title,format='%.0e')
    plt.xlabel('X position (m)')
    plt.ylabel('Y position (m)')
