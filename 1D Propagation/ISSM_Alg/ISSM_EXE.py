import numpy as np
from ISSM_Alg.SMM_conversion import convert_1D_SMM,convert_k_SMM
from ISSM_Alg.Matrices_class_set_SMM import scattering_matrix_set

'''
def execute_ISSM(ray_data,n_real,n_imag,antenna,freq):
    ray_1D,ret_value = convert_1D_SMM(ray_data,antenna)
    sms = scattering_matrix_set()
    if ret_value==0:
        ray_k_values = convert_k_SMM(n_real,n_imag,freq)
        coeff_ray = sms.Total_coefficients(ray_k_values,ray_1D)
        #print(coeff_ray)
        A_mag = np.abs((coeff_ray[0]))
        A_phase = np.arctan(np.imag(coeff_ray[0])/np.real(coeff_ray[0]))
        D_mag = np.abs((coeff_ray[1])) #* self.n_real[-1]/self.n_real[0]
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
    """
    Execute the Integrated Scattering Signal Model (ISSM) for a layered plasma structure.

    This function computes the transmission, reflection, and intermediate layer
    coefficients for an RF signal propagating through a medium with complex 
    refractive indices.

    Parameters:
    length : float
        Total length (thickness) of the stratified medium (in meters).
    n_real : list or ndarray of floats
        Real part of the refractive index in each layer.
    n_imag : list or ndarray of floats
        Imaginary part of the refractive index in each layer.
    freq : float
        Frequency of the incident wave (in Hz).

    Returns:
    A_mag : float
        Magnitude of the forward-propagating wave at the input (reflection coefficient magnitude).
    D_mag : float
        Magnitude of the transmitted wave at the output (transmission coefficient magnitude).
    A_phase : float
        Phase angle (in radians) of the forward wave at the input.
    D_phase : float
        Phase angle (in radians) of the transmitted wave at the output.
    Coeff_B : list of complex
        Partial reflection coefficients at intermediate interfaces.
    Coeff_C : list of complex
        Partial transmission coefficients at intermediate interfaces.
    B_phase : list of float
        Phase angles (radians) for intermediate reflected waves.
    C_phase : list of float
        Phase angles (radians) for intermediate transmitted waves.
    
    Notes:
    - The function assumes the existence of a `scattering_matrix_set` class
      that provides `Total_coefficients` and `partial_coefficients` methods.
    - The wave vector components are computed from complex refractive indices
      using the `convert_k_SMM` function.
    """
    sms = scattering_matrix_set()
    ray_k_values = convert_k_SMM(n_real,n_imag,freq)
    coeff_ray = sms.Total_coefficients(ray_k_values,length)
    Coeff_B, Coeff_C,B_phase,C_phase = sms.partial_coefficients(ray_k_values,length,coeff_ray[0],coeff_ray[1])
    #print(coeff_ray)
    A_mag = np.abs((coeff_ray[0]))**2
    A_phase = np.arctan(np.imag(coeff_ray[0])/np.real(coeff_ray[0]))
    D_mag = np.abs((coeff_ray[1]))**2 #* n_real[-1]/n_real[0]
    D_phase = np.arctan(np.imag(coeff_ray[1])/np.real(coeff_ray[1]))
    #Abs = 1-A_mag**2 - D_mag**2
    return A_mag, D_mag, A_phase, D_phase,Coeff_B,Coeff_C,B_phase,C_phase