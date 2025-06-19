import numpy as np
from dataclasses import dataclass

@dataclass
class input_data:
    #####################################################
    # Physical Parameters ###############################
    #####################################################
    '''
    Not all parameters here will be used. The used parameters depends on the selected data type ran in the main propagation file. 
    '''    
    
    #Juan wedge: 
    cfd_data_file = "M15_30km.csv"
    X_range = [-0.052810,-0.005040]
    X_input = -0.025192#-0.01061#-0.035
    Y_input = 0.0095622#0.011372#0.008364
    Y_range = [0.006206,0.0120672]
    #coord points: (-0.035,0.0131), (-0.02,0.0301), (-0.05,0.00968)
    '''
    cfd_data_file = "Heating_Up.csv" #name of CFD file with final data
    X_range = [0.14679 ,0.50079] #[-0.06,-0.04742] # Linear bounds of field data (X-component)  
    X_input = 0.25037 # Point of propagation
    Y_input = 0.06555# Heating_up 0.06675# Muadib 0.06555 
    Y_range = [0.047758,0.11398]  #[0.0055325,0.0070775] # Linear bounds of field data (Y-component)
    '''
    height = 0.018 # Length of propagation Muadib 0.1 Juan 0.018 Heating_up 0.084
    Frequency = 14e9#np.linspace(11e9,14e9,100) #frequency of propagating rays (Hz)
    frequency_range = [5e9, 8e9, 11e9, 14e9]
    B_field = [0.0,0.0,1.0] #magnetic field vector
    B_range = np.linspace(0.0,1.0,100) # Range of Magnetic field magnitudes
    #####################################################
    # Computation #######################################
    #####################################################
    
    run_type = 'FComp'
    