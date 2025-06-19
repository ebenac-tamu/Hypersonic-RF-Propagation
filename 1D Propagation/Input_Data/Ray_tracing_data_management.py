import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
import csv
from matplotlib.colors import LogNorm
from scipy import interpolate
from Input_Data.Hypersonic_Collision_Frequency import calculate_collisions
from globalz import global_variables 
from matplotlib.colors import SymLogNorm
#from Full_Execute import csv_file
#from INPUTS import input_data as i_d
#import functools
glbl = global_variables()

class import_data_refractive_index:
    #@functools.cached_property
    def __init__(self,input_data):
        print('Initializing data...')
        """
        Initialize and load CFD data from a CSV file, compute electron-neutral collision frequency,
        and create interpolators for spatial access to key plasma properties.

        Parameters:
            input_data (str): Path to the CSV file containing CFD simulation data.
        
        Stores:
            - Electron density, vibrational temperature, air density.
            - Densities of various neutral species (O, O2, N, N2, NO).
            - Total neutral density.
            - Electron-neutral collision frequency.
            - Interpolators for electron density, neutral refraction, and collision frequency.
        """
        
        oN = 1.1
        oN2 = 1.71
        oO = 0.802
        oO2 = 1.562
        oNO = 1.698
        ev = np.array([0.001,0.0015,0.0018,0.002,0.0025,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.012,0.015,0.018,0.02,0.025,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.15,0.18,0.2,0.25,0.3,0.4,0.5])
        cross_section = np.array([1.357,1.426,1.464,1.49,1.55,1.62,1.718,1.81,1.908,2.0,2.062,2.131,2.19,2.342,2.55,2.729,2.85,3.12,3.4,3.85,4.33,4.72,5.1,5.41,5.69,5.95,6.45,7.1,7.59,7.9,8.5,9.0,9.7,10.16])*1e-20
        cs = interpolate.CubicSpline(ev,cross_section)
        print('Extracting CFD Data...')
        data_set = open(input_data,'r')
        file = csv.DictReader(data_set)
        self.x_pos = []
        self.y_pos = []

        n_O = []
        n_O2 = [] 
        n_N = []
        n_N2 = []
        n_NO = []
        N_O = []
        N_O2 = [] 
        N_N = []
        self.N_N2 = []
        N_NO = []
        self.ele = []
        Temp_vib = []
        rho_air = []
        #self.nu = []
        for row in file:
            #print(row)
            self.x_pos.append(float(row['x']))
            self.y_pos.append(float(row['y']))
            n_O.append(float(row['refra_O=']))
            n_O2.append(float(row['refra_O2/RO_O2='])) 
            n_N.append(float(row['REFRA_N='])) 
            n_N2.append(float(row['REFRA_N2=']))
            n_NO.append(float(row['REFRA_NO=']))
            self.ele.append(float(row['nD_e-']))
            Temp_vib.append(float(row['Tv']))
            rho_air.append(float(row['rho']))
            N_N.append(float(row['nD_N']))
            self.N_N2.append(float(row['nD_N2']))
            N_O.append(float(row['nD_O']))
            N_O2.append(float(row['nD_O2']))
            N_NO.append(float(row['nD_NO']))
            #self.nu.append(1e10)#float(row['nu']))
        print('CFD Data Extracted')
        print('Calculating Collision Frequency...')

        self.N_neutral = np.zeros_like(rho_air)  
        self.n_neutral = np.zeros_like(rho_air)  
        self.nu = np.zeros_like(rho_air)    
        #self.nu = np.zeros_like(n_neutral)
        for i in range(len(rho_air)):
            self.n_neutral[i] = n_N[i]+n_N2[i]+n_O[i]+n_O2[i]+n_NO[i]#2e-30 * np.pi *(N_N[i]*oN+self.N_N2[i]*oN2+N_O[i]*oO+N_O2[i]*oO2+N_NO[i]*oNO)
            self.N_neutral[i] = N_N[i]+self.N_N2[i]+N_O[i]+N_O2[i]+N_NO[i]
            self.nu[i] = 5e9#cs(Temp_vib[i]/11600) * self.N_neutral[i] * np.sqrt(Temp_vib[i]*glbl.Kb/glbl.m_e)
            self.ele[i] = self.ele[i]
        #n_neutral = np.nan_to_num(n_neutral)
        '''
        with open("nu.txt", "w") as f:
            f.write("(\n")
            for value in self.nu:
                f.write(f"{value}\n")
            f.write(")")
        '''
        print('Collision Frequency Calculated')
        #print(len(n_neutral))
        #print(len(nu))
        print('Interpolating Data...')
        self.refra_neutral = interpolate.CloughTocher2DInterpolator(list(zip(self.x_pos,self.y_pos)),self.n_neutral)
        self.collision_frequency = interpolate.CloughTocher2DInterpolator(list(zip(self.x_pos,self.y_pos)),self.nu)
        self.electron_density = interpolate.CloughTocher2DInterpolator(list(zip(self.x_pos,self.y_pos)),self.ele)       
        print('Data Interpolated')
        #self.cache = {}
        print('Data initialized')
    
    def input_variable_parameters(self,frequency,B_field,height=0,m=0,b=-1000):
        """
        Define key simulation input parameters like frequency, magnetic field, and vertical mesh.
        
        Parameters:
            frequency (float): Frequency of the EM wave (Hz).
            B_field (np.array): Magnetic field vector (Tesla).
            height (float, optional): Vertical domain size to generate height profile (default 0).
            m (float, optional): Slope of the ray trajectory (default 0).
            b (float, optional): Intercept of the ray trajectory (default -1000).
        """
        self.m = m
        self.b = b
        self.freq = frequency
        self.B_mag = np.sqrt(B_field[0]**2+B_field[1]**2+B_field[2]**2)
        self.B_unit = B_field/self.B_mag
        if height != 0:
            self.h = np.linspace(00.00,height,1000)
    
    def refractive_index_set(self,x_val,y_val,direct):
        """
        Compute total refractive index and its spatial gradient at a specific (x, y) location.

        Parameters:
            x_val (float): x-coordinate.
            y_val (float): y-coordinate.
            direct (not used): Placeholder for direction (not implemented).

        Returns:
            tuple: (n, dndx, dndy)
                n (float): Total refractive index at the point.
                dndx (float): Gradient of refractive index in x-direction.
                dndy (float): Gradient of refractive index in y-direction.
        """
        #distance = (y_val - (self.m * x_val + self.b))
        distance = 1
        if np.real(distance) >= 0:
            x_val = np.real(x_val)
            y_val = np.real(y_val)
            step = 0.0001
            B_dir = tuple(self.B_unit)
            direction = (direct[0],direct[1],0)
            #angle = glbl.vector_angle(B_dir,direction)
            angle = np.pi/2
            
            n_neutral_0 = self.refra_neutral(x_val,y_val)
            n_neutral_mx1 = self.refra_neutral(x_val-step,y_val)
            n_neutral_py1 = self.refra_neutral(x_val,y_val+step)
            
            Ne_0 = self.electron_density(x_val,y_val)
            Ne_mx1 = self.electron_density(x_val-step,y_val)
            Ne_py1 = self.electron_density(x_val,y_val+step)

            nu_0 = self.collision_frequency(x_val,y_val)
            nu_mx1 = self.collision_frequency(x_val-step,y_val)
            nu_py1 = self.collision_frequency(x_val,y_val+step)
            
            n_plasma_0 = np.sqrt(glbl.Appleton_Hartree(Ne_0,self.B_mag,self.freq,nu_0,angle,polar='lhc'))-1
            #n_plasma_0 = np.nan_to_num(n_plasma_0)
            n_plasma_mx1 = np.sqrt(glbl.Appleton_Hartree(Ne_mx1,self.B_mag,self.freq,nu_mx1,angle,polar='lhc'))-1
            n_plasma_py1 = np.sqrt(glbl.Appleton_Hartree(Ne_py1,self.B_mag,self.freq,nu_py1,angle,polar='lhc'))-1
            #print(n_plasma_0)
            n_true_0 = n_plasma_0+1+n_neutral_0
            #print(n_true_0)
            n_true_mx1 = n_plasma_mx1+1 +n_neutral_mx1
            n_true_py1 = n_plasma_py1+1 +n_neutral_py1
            
            dndx = (n_true_mx1-n_true_0)/step
            dndy = (n_true_py1-n_true_0)/step
        elif np.real(distance) <= 0:
            norm = np.sqrt(self.m**2 + 1)
            n_true_0=100
            dndx = -100 * (self.m / norm)  
            dndy = 100 * (1 / norm)  
        return n_true_0, dndx, dndy
    
    def extract_data_for_plotting(self):
        """
        Retrieve CFD spatial coordinates and electron density for external plotting.

        Returns:
            tuple: (x_pos, y_pos, electron_density)
        """
        return self.x_pos,self.y_pos,self.ele
    
    def extract_SMM_inputs(self,x_value,y_initial,ij =0):
        """
        Extract refractive index profiles (real and imaginary parts) along a vertical column
        for different collision frequency scaling cases.

        Parameters:
            x_value (float): x-position.
            y_initial (float): y starting point (bottom of profile).
            ij (int, optional): Mode selector:
                - 0: Nominal collision frequency.
                - 1: 5x collision frequency.
                - 2: 0.2x collision frequency.

        Returns:
            tuple: (n_real, n_imag, h)
                n_real (list): Real part of total refractive index.
                n_imag (list): Imaginary part of total refractive index (absorption).
                h (np.array): Height vector used for profiling.
        """
        n_real = []
        n_imag = []
        #y_initial = 0.08798313
        for i in range(len(self.h)):
            #print(self.collision_frequency(x_value,self.h[i]+y_initial))
            if ij ==0:
                n_val = np.sqrt(glbl.AH_parallel(self.electron_density(x_value,self.h[i]+y_initial),self.B_mag,self.freq,self.collision_frequency(x_value,self.h[i]+y_initial)))
            elif ij ==1:
                n_val = np.sqrt(glbl.AH_parallel(self.electron_density(x_value,self.h[i]+y_initial),self.B_mag,self.freq,self.collision_frequency(x_value,self.h[i]+y_initial)*5))
            elif ij ==2:
                n_val = np.sqrt(glbl.AH_parallel(self.electron_density(x_value,self.h[i]+y_initial),self.B_mag,self.freq,self.collision_frequency(x_value,self.h[i]+y_initial)/5))
            n_imag.append(-1*np.imag(n_val))
            n_real.append(np.real(n_val)+self.refra_neutral(x_value,self.h[i]+y_initial))
        return n_real,n_imag,self.h
    
    def Ratio_data(self,x_initial,y_initial):
        """
        Compute plasma frequency to wave frequency ratio (X) and
        normalized collision frequency (Y) across height.

        Parameters:
            x_initial (float): x-coordinate of evaluation.
            y_initial (float): y-coordinate of the bottom of the profile.

        Returns:
            tuple: (X, Y, h)
                X (list): Plasma frequency / wave frequency ratio.
                Y (list): Collision frequency / wave frequency ratio.
                h (np.array): Height vector.
        """
        X = []
        Y = []
        for i in range(len(self.h)):
            X.append(glbl.plasma_frequency(self.electron_density(x_initial,y_initial+self.h[i]))/(self.freq*2*glbl.pi))
            Y.append(self.collision_frequency(x_initial,y_initial+self.h[i])/(self.freq*2*glbl.pi))
        return X,Y,self.h

    def plot_cfd_data(self):
        """
        Visualize different CFD-derived plasma parameters:
            - Collision frequency.
            - Refractive index (real part minus 1).
            - Refractive index (imaginary part).
            - Plasma-to-collision frequency ratio.
        
        Note:
            Uses `matplotlib.pyplot` to generate scatter plots of spatial data.
        """
        '''
        plt.scatter(self.x_pos, self.y_pos, c=self.N_neutral,s=60, cmap='jet', alpha=0.5, norm=LogNorm(vmin=5e23,vmax=2e24))  # 'c' argument specifies the color
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='Neutral Density ($m^{-3}$)',format='%.0e')  # Add color bar
        plt.show()
        plt.scatter(self.x_pos, self.y_pos, c=self.n_neutral,s=20, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e-6,vmax=1e-3))  # 'c' argument specifies the color
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='n_neutral',format='%.0e')  # Add color bar
        plt.show()
        plt.scatter(self.x_pos, self.y_pos, c=self.ele,s=20, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e16,vmax=1e19))  # 'c' argument specifies the color
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='Electron Density ($m^{-3}$)',format='%.0e')  # Add color bar
        plt.show()
        '''
        plt.scatter(self.x_pos, self.y_pos, c=self.nu,s=20, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e7,vmax=1e11))  # 'c' argument specifies the color
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='Collision Frequency ($s^{-1}$)',format='%.0e')  # Add color bar
        plt.show()
        
        
        n_real = []
        n_imag = []
        #X = []
        ratio = []
        for i in range(len(self.ele)):
            n_val = np.sqrt(glbl.AH_parallel(self.ele[i],self.B_mag,self.freq,self.nu[i]))+self.n_neutral[i]
            n_imag.append(-1*np.imag(n_val))
            #X.append(glbl.plasma_frequency(self.ele[i])/(self.freq*2*np.pi))
            ratio.append(glbl.plasma_frequency(self.electron_density(self.x_pos[i],self.y_pos[i]))**2/(self.nu[i]*self.freq))
            #print(np.real(n_val))
            n_real.append(np.real(n_val)-1)
        plt.scatter(self.x_pos, self.y_pos, c=ratio,s=3, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e-1,vmax=1e3))  # 'c' argument specifies the color
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='Y',format='%.0e')  # Add color bar
        plt.show()
        #plt.scatter(self.x_pos, self.y_pos, c=n_real,s=40, cmap='viridis', alpha=0.5)  # 'c' argument specifies the color
        #plt.xlabel('X-axis')
        #plt.ylabel('Y-axis')
        #plt.colorbar(label='n_real',format='%.0e')  # Add color bar
        #plt.show()
        plt.scatter(self.x_pos, self.y_pos, c=n_imag,s=3, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e-5,vmax=1e-1))  # 'c' argument specifies the color
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='Im(n)',format='%.0e')  # Add color bar
        plt.show()
        plt.scatter(self.x_pos, self.y_pos, c=n_real,s=3, cmap='jet', alpha=0.5,norm=SymLogNorm(vmin=-1e0,vmax=1e0,linthresh=1e-1)) # real scatter plot
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='Re(n-1)',format='%.0e')  # Add color bar
        plt.show()

    def plot_vertical_CS(self,X_point,Y_point):
        plotting_data = []
        for i in range(len(self.height)):
            plotting_data.append(self.electron_density())

        return 0 
    
    def ray_plot_background(self): 
        """
        Plot electron density field as a background (e.g., for ray tracing overlay).
        
        Returns:
            None (plots directly using matplotlib).
        """
        plt.scatter(self.x_pos, self.y_pos, c=self.ele,s=20, cmap='jet', alpha=0.5, norm=LogNorm(vmin=1e16,vmax=1e19))  # 'c' argument specifies the color
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.colorbar(label='Electron Density ($m^{-3}$)',format='%.0e')  # Add color bar
    
    

        
#idri = import_data_refractive_index(i_d.cfd_data_file)