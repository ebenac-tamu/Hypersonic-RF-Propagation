#from Input_Data.Ray_tracing_data_management import import_data_refractive_index
from Input_Data.Ray_Extraction_old import Extract_Ray_Data
from Ray_Tracing_Alg.eikonal import ray_tracing
import numpy as np

def execute_RTA(Initial_X,Initial_Y,artificial_antenna,rf,angle_range = [np.pi/3,np.pi*7/8],rays=50,refinement=True,step = 0.00005):
    Initial_angle_range = np.linspace(angle_range[0],angle_range[1],rays)
    aa_x1,aa_y1 = artificial_antenna[0]
    aa_x2,aa_y2 = artificial_antenna[1]
    x_aa_max = max([aa_x1,aa_x2])
    x_aa_min = min([aa_x1,aa_x2])
    y_aa_min = Initial_Y+0.002
    x_max = max([aa_x1,aa_x2])+0.0005
    x_min = min([aa_x1,aa_x2])-0.0005
    y_max = max([aa_y1,aa_y2])+0.0000001
    y_min = min([aa_y1,aa_y2])-0.07
    
    X_lim = (x_min,x_max)
    Y_lim = (y_min,y_max)
    rte = ray_tracing(X_lim,Y_lim,rf)
    #x_data_set = []
    #y_data_set = []
    ray_data = []
    n_real_set = []
    n_imag_set = []
    if refinement==True:
        print('Locating Refinement Region...')
        theta_new = []
        for theta in Initial_angle_range:
            x_fin,y_fin = rte.eikonal_fast(Initial_X,Initial_Y,theta)
            #print(aa_x1)
            #print(x_fin)
            #print(aa_x2)
            #print()
            if x_aa_min<x_fin<x_aa_max and y_fin>y_aa_min:
                #print(theta)
                theta_new.append(theta)
        #print(theta_new)
        theta_new_min = min(theta_new)
        theta_new_max = max(theta_new)
        Initial_angle_range = np.linspace(theta_new_min,theta_new_max,rays)
        print('Angular region refined to range from '+str(theta_new_min*180/np.pi)+' to '+str(theta_new_max*180/np.pi))
    counter = 0
    for theta in Initial_angle_range:
        counter +=1
        print('Calculating ray number '+str(counter))
        x_dat,y_dat = rte.eikonal_fast(Initial_X,Initial_Y,theta,data_rec=True)
        #print(ray)
        #print(x_data)
        #x_data_set.append(x_data)
        #y_data_set.append(y_data)
        ray = []
        for i in range(len(x_dat)):
            ray.append(tuple((x_dat[i],y_dat[i])))
        ray_data.append(ray)
        n_real=[]
        n_imag=[]
        for i in range(len(x_dat)):
            n_comp,nul,nul1 = rf.refractive_index_set(x_dat[i],y_dat[i],[1,0])
            n_real.append(np.real(n_comp))
            n_imag.append(np.imag(n_comp))
        n_real_set.append(n_real)
        n_imag_set.append(n_imag)
    
    #erd = Extract_Ray_Data(rays,x_data_set,y_data_set)
    #erd.ray_data_extract()
    print('Ray data calculated')
    return ray_data,n_real_set,n_imag_set

    
    
    
    
    