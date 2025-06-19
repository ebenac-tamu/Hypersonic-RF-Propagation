import numpy as np
from globalz import global_variables
glbl = global_variables()

def convert_1D_SMM(ray_list,antenna):
    """
    Converts a 2D ray path into a 1D spatial mapping along the propagation direction
    and checks if the ray exits within antenna bounds.

    Parameters:
    ray_list : list of tuples
        List of (x, y) coordinate tuples representing the ray path.
    antenna : tuple of tuples
        Contains two points (each a tuple of x, y) representing the antenna aperture ends.

    Returns:
    Ray_R : list of floats
        Cumulative distance from the ray origin at each point along the path.
        This effectively maps the ray path into a 1D length vector.
    ret : int
        A binary flag indicating whether the ray exits within the antenna bounds:
        - 0: Ray terminates within antenna x-range and above x-axis (valid transmission).
        - 1: Ray misses the antenna aperture or terminates below x-axis (invalid transmission).
    """
    x_aa_1,y_aa_1 = antenna[0]
    x_aa_2,y_aa_2 = antenna[1]
    x_min = min([x_aa_1,x_aa_2])
    x_max = max([x_aa_1,x_aa_2])
    ray_primary = ray_list
    # print(ray_primary)
    Ray_R = [0]
    length = 0
    for i in range(0,len(ray_primary)-1):
        #print(ray_primary[i])
        P1_x1,P1_y1 = ray_primary[i]
        P1_xp1,P1_yp1 = ray_primary[i+1]
        length += np.sqrt((P1_x1-P1_xp1)**2+(P1_y1-P1_yp1)**2)
        Ray_R.append(length)
    x_fin, y_fin =ray_primary[-1]
    print(len(Ray_R))
    print(Ray_R)
    if x_min<x_fin<x_max and y_fin>0.0:
        
        ret = 0
        #print(ret)
    else:
        ret = 1
    return Ray_R, ret

def convert_k_SMM(n_real,n_imag,freq):
    """
    Converts complex refractive index profiles into spatially-resolved wave numbers 
    (propagation constants) for scattering matrix method (SMM) calculations.

    Parameters:
    n_real : list or array-like of float
        Real parts of the refractive index at discrete spatial locations.

    n_imag : list or array-like of float
        Imaginary parts of the refractive index at the same locations.
        This typically accounts for attenuation/absorption in the medium.

    freq : float
        Frequency of the electromagnetic wave (Hz).

    Returns:
    k_list : list of complex
        List of complex wave numbers (propagation constants) computed from the 
        complex refractive indices. Each value corresponds to a spatial segment.
    """
    k_list = []
    for i in range(len(n_real)):
        n = n_real[i]-1j*n_imag[i]
        k_list.append(glbl.propagation_constant(freq,n))
    return k_list