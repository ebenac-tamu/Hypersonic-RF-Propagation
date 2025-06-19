import numpy as np
from ISSM_Alg.SMM_conversion import convert_1D_SMM,convert_k_SMM
from ISSM_Alg.Matrices_class_set_SMM import scattering_matrix_set

def execute_ISSM(ray_data,n_real,n_imag,antenna,freq):
    ray_1D,ret_value = convert_1D_SMM(ray_data,antenna)
    sms = scattering_matrix_set()
    if ret_value==0:
        ray_k_values = convert_k_SMM(n_real,n_imag,freq)
        coeff_ray = sms.Total_coefficients(ray_k_values,ray_1D)
        #print(coeff_ray)
        A_mag = np.abs((coeff_ray[0]))**2
        A_phase = np.arctan(np.imag(coeff_ray[0])/np.real(coeff_ray[0]))
        D_mag = np.abs((coeff_ray[1]))**2 #* self.n_real[-1]/self.n_real[0]
        D_phase = np.arctan(np.imag(coeff_ray[1])/np.real(coeff_ray[1]))
        #Abs = 1-A_mag**2 - D_mag**2
    elif ret_value==1:
        A_mag = 0
        A_phase = 0
        D_mag = 0
        D_phase = 0
        Abs = 0
        #Coeff_B = 0
        #Coeff_C = 0 
        #B_phase = 0
        #C_phase = 0
    return A_mag, D_mag, A_phase, D_phase#,Coeff_B,Coeff_C,B_phase,C_phase
'''
def execute_ISSM(length,n_real,n_imag,freq):
    sms = scattering_matrix_set()
    ray_k_values = convert_k_SMM(n_real,n_imag,freq)
    coeff_ray = sms.Total_coefficients(ray_k_values,length)
    #print(coeff_ray)
    A_mag = np.abs((coeff_ray[0]))**2
    A_phase = np.arctan(np.imag(coeff_ray[0])/np.real(coeff_ray[0]))
    D_mag = np.abs((coeff_ray[1]))**2 #* n_real[-1]/n_real[0]
    D_phase = np.arctan(np.imag(coeff_ray[1])/np.real(coeff_ray[1]))
    #Abs = 1-A_mag**2 - D_mag**2
    return A_mag, D_mag, A_phase, D_phase#,Coeff_B,Coeff_C,B_phase,C_phase
'''