import numpy as np
from dataclasses import dataclass

@dataclass
class input_data:
    #####################################################
    # Physical Parameters ###############################
    #####################################################
    #X_initals = [0.152,0.2038,0.2515,0.3,0.35,0.4,0.4527,0.502,0.5466,0.6,0.6456,0.694]
    #Y_initals = [0.1505,0.1585,0.1657,0.1732,0.181,0.1886,0.1966,0.2042,0.2112,0.21955,0.2265,0.2339]
    #Surface Points List: (0,0) (0.001095,0.00287) (0.004956,0.0050812) (0.01011,0.005733) (0.02002,0.007006) (0.03022,0.008291) (0.04013,0.009539) (0.05007,0.010815) (0.06003,0.012075) (0.064904,0.012714)
    cfd_data_file = "muadib.csv" #name of CFD file with final data
    X_initial = 0.25037 #initial X-position of rays
    Y_initial = 0.06555 #initial Y-position of rays
    Frequency = 1.3e9 #frequency of propagating rays (Hz)
    B_field = [0.0,0.0,0.0] #magnetic field vector
    Artificial_Antenna = [(X_initial+0.01 ,Y_initial+0.08),(X_initial-0.01,Y_initial+0.08)] #location of artificial antenna in 2D [(x1,y1),(x2,y2)]
    
    #####################################################
    # Computational Parameters ##########################
    #####################################################
    angle_range = [84*np.pi/180,94*np.pi/180]
    rays = 50
    refinement = False
    raytracing_step_size = 0.000005
    

