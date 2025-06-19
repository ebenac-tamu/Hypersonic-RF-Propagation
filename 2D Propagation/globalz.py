import numpy as np

class global_variables:
    def __init__(self):
        self.c = 2.99792458e8 # Speed of Light (m/s)
        self.mu_0 = np.pi*4e-7 # Vacuum Magnetic Permiability (N*A^-2)
        self.m_e = 9.11e-31 # Electron mass (kg)
        self.q = 1.6e-19 # Elementary Charge (C)
        self.esp_0 = 8.85e-12 # Vacuum Electric Permativitty (F/m)
        self.h = 6.626e-34 # Plank Constant (J*s)
        self.Kb = 1.38e-23 # Boltzmann Constant (J/K)
        self.Na = 6.022e23 # Avagadros Number (particles/mole)
        self.R = 8.3145 # Universal Gas Constant (J/(mol*K))
        self.F = 96485.3342 # Faradays Constant (C/mol)
        self.G = 6.674e-11 # Gravitational Constant (m^3/(kg*s^2))
        self.pi = 3.1415926535 # Pi
        self.e = 2.7182818284 # Eulers Number
        self.sigma = 5.67e-8 # Stefan-Boltzmann Constant (W/(m^2 K^4))

    ######################################################################
    ######################################################################
    # OPTICS #############################################################
    ######################################################################
    ######################################################################  

    def propagation_constant(self,f_0,refra):
        '''
        f_0 - scalar - electromagnetic frequency (Hz)
        refra - scalar - refractive index
        '''
        W_0 = f_0*2*np.pi
        return W_0/self.c * refra
    
    def critical_angle(self,n1,n2):
        # n1 - scalar - initial refractive index
        # n1 - scalar - final refractive index
        return np.arcsin(n2/n1)
    
    def beer_lambert(self,f_0,refra,del_x):
        '''
        f_0 - scalar - electromagnetic frequency (Hz)
        refra - scalar - refractive index
        del_x - scalar - propagation displacement (m)
        '''
        w_0 = f_0*(2*np.pi)
        Xi = -1*np.imag(refra)
        return (np.exp(-1*Xi*(w_0/self.c)*del_x))**2

    ######################################################################
    ######################################################################
    # PLASMA PROPERTIES ##################################################
    ######################################################################
    ######################################################################

    def plasma_frequency(self,Ne):
        # Ne - scalar - electron numbder density (m^-3)
        return np.sqrt((Ne*self.q**2)/(self.m_e*self.esp_0))
    
    def critical_electron_density(self,f_0):
        # f_0 - scalar - electromagnetic frequency (Hz)
        W_0 = f_0*2*np.pi
        return W_0**2 * self.m_e*self.esp_0/(self.q**2)

    def gyrofrequency(self,B):
        # B - scalar - magnetic field magnitude (T)
        return self.q*B/(self.m_e)
    
    def electron_thermal_velocity(self,Te):
        # most probable thermal velocity
        # Te - scalar - electron temperature (K)
        return np.sqrt(8*self.Kb*Te/(self.m_e*np.pi))
    
    def debye_length(self,Te,Ne):
        # Te - scalar - eletron temperature (K)
        # Ne - scalar - electron numbder density (m^-3)
        return np.sqrt(self.esp_0*self.Kb*Te/(Ne*self.q**2))
    
    def ExB_drift(self,E_field,B_field):
        # E_field - vector - Electric field (V/m)
        # B_field - vector - Magnetic field (T)
        return np.cross(E_field,B_field)/np.abs(B_field)
    
    def Hall_Parameter(self,B_mag,nu):
        # B_mag - scalar - Magnetic field strength (T)
        # nu - scalar - collision frequency (1/s)
        return self.q*B_mag/(self.m_e*nu)
    
    
    
    ######################################################################
    ######################################################################
    # PLASMA KINETICS ####################################################
    ######################################################################
    ######################################################################

    def coulomb_collision_frequency(self,Ne,Te):
        # Ne - scalar - electron numbder density (m^-3)
        # Te - scalar - eletron temperature (K)
        coul_log = 0.5*np.log(1.0+0.25*self.esp_0*self.Kb*Te/(np.pi*self.q**2*Ne*(self.q**4/(16*np.pi**2*self.esp_0**2*self.Kb**2*Te**2)+self.h**2/(2*np.pi*self.Kb*self.m_e*Te))))
        return 3.633e-6*coul_log*Ne*Te**(-1.5)
    
    def electron_mobility(self,nu):
        # nu - scalar - collision frequency (1/s)
        return self.q/(self.m_e*nu)
    
    def electron_diffusion_coefficient(self,Temp,nu):
        # Temp - scalar - gas temperature (K)
        # nu - scalar - collision frequency (1/s)
        return self.Kb*Temp/(self.m_e*nu)
    
    def bohm_diffusion_coefficient(self,B_mag,Te):
        # Te - scalar - electron temperature (K)
        # B_mag - scalar - magnitude of magnetic field strength (T)
        return self.Kb*Te/(16*self.q*B_mag)

    ######################################################################
    ######################################################################
    # PLASMA OPTICS ######################################################
    ######################################################################
    ######################################################################

    def AH_parallel(self,Ne,B_mag,f_0,nu):
        # Ne - scalar - electron numbder density (m^-3)
        # B_mag - scalar - magnetic field magnitude (T)
        # f_0 - scalar - electromagnetic frequency (Hz)
        # nu - scalar - collision frequency (s^-1)
        omega_p = self.plasma_frequency(Ne)
        omega_c = self.gyrofrequency(B_mag)
        frequency = f_0*2*np.pi
        n2 = 1-(omega_p/frequency)**2/(1-1j*(nu/frequency)-(omega_c/frequency)**2/(1-1j*(nu/frequency)-(omega_p/frequency)**2))
        return n2

    def AH_perpendicular(self,Ne,B_mag,f_0,nu):
        # Ne - scalar - electron numbder density (m^-3)
        # B_mag - scalar - magnetic field magnitude (T)
        # f_0 - scalar - electromagnetic frequency (Hz)
        # nu - scalar - collision frequency (s^-1)
        omega_p = self.plasma_frequency(Ne)
        omega_c = self.gyrofrequency(B_mag)
        frequency = f_0*2*np.pi
        n2 = 1-(omega_p/frequency)**2/(1-1j*(nu/frequency)-(omega_c/frequency))
        return n2

    def Appleton_Hartree(self,Ne,B_mag,f_0,nu,angle,polar = 'rhc'):
        # Ne - scalar - electron numbder density (m^-3)
        # B_mag - scalar - magnetic field magnitude (T)
        # f_0 - scalar - electromagnetic frequency (Hz)
        # nu - scalar - collision frequency (s^-1)
        # angle - scalar - angle between magnetic field and propagation constant (radians)
        # polar - string - left/right hand circular polarization
        omega_p = self.plasma_frequency(Ne)
        omega_c = self.gyrofrequency(B_mag)
        frequency = f_0*2*np.pi
        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)
        #X = (omega_p/frequency)**2
        #Y = omega_c/frequency
        #Z = nu/frequency
        #val1 = 1-1j*Z-X
        #val2 = 0.5*Y**2*sin_theta**2
        denom_1 = 1-1j*(nu/frequency)-(0.5*(omega_c/frequency)**2 * sin_theta**2)/(1-1j*(nu/frequency)-((omega_p/frequency)**2))
        denom_2 = np.sqrt((0.25*(omega_c/frequency)**4 * sin_theta**4)+(omega_c/frequency)**2 * cos_theta**2 * (1-1j*(nu/frequency)-((omega_p/frequency)**2))**2)/(1-1j*(nu/frequency)-((omega_p/frequency)**2))
        if polar == 'rhc':
            n2 = 1-((omega_p/frequency)**2)/((1-1j*(nu/frequency)-(0.5*(omega_c/frequency)**2 * sin_theta**2)/(1-1j*(nu/frequency)-((omega_p/frequency)**2)))-(np.sqrt((0.25*(omega_c/frequency)**4 * sin_theta**4)+(omega_c/frequency)**2 * cos_theta**2 * (1-1j*(nu/frequency)-((omega_p/frequency)**2))**2)/(1-1j*(nu/frequency)-((omega_p/frequency)**2))))
            #n2 = (B-np.sqrt(B**2-4*A*C))/(2*A)
        elif polar == 'lhc':
            n2 =1-((omega_p/frequency)**2)/(denom_1+denom_2)
            #n2 = (B+np.sqrt(B**2-4*A*C))/(2*A)
        
        return n2
    
    def cold_cutoff(self,Ne,B_mag,f_0):
        # Ne - scalar - electron numbder density (m^-3)
        # B_mag - scalar - magnetic field magnitude (T)
        # f_0 - scalar - electromagnetic frequency (Hz)
        W_p = self.plasma_frequency(Ne)
        W_c = self.gyrofrequency(B_mag)
        W_0 = f_0*2*np.pi
        if W_p**2-(W_c+W_0)<=0.00000001 or W_p**2-(W_0-W_c)<=0.00000001: 
            return True
        else:
            return False
        
    def cold_resonance(self,Ne,B_mag,f_0):
        # Ne - scalar - electron numbder density (m^-3)
        # B_mag - scalar - magnetic field magnitude (T)
        # f_0 - scalar - electromagnetic frequency (Hz)
        W_p = self.plasma_frequency(Ne)
        W_c = self.gyrofrequency(B_mag)
        W_0 = f_0*2*np.pi
        return W_p**2-(W_0**2-W_c**2)<=0.00000001

    def air_to_plasma_transmission(self,f_0,n_plasma):
        # f_0 - scalar - electromagnetic frequency (Hz)
        # n_plasma - scalar - refractive index of plasma
        K_0 = self.propagation_constant(f_0,1.0)
        K_1 = self.propagation_constant(f_0,n_plasma)
        return 2*K_0/(K_0+K_1)
    
    def air_to_plasma_reflection(self,f_0,n_plasma):
        # f_0 - scalar - electromagnetic frequency (Hz)
        # n_plasma - scalar - refractive index of plasma
        K_0 = self.propagation_constant(f_0,1.0)
        K_1 = self.propagation_constant(f_0,n_plasma)
        return (K_0-K_1)/(K_0+K_1)

    ######################################################################
    ######################################################################
    # KINETICS ###########################################################
    ######################################################################
    ######################################################################

    def N2_cross_section_low_ev(self,Te):
        # Te - scalar - electron temperature below 4ev (K)
        Te = Te/11600
        ev_data = np.array([0.0015,0.0018,0.002,0.0025,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.012,0.015,0.018,0.02,0.025,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.15,0.18,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,1.92,1.98,2.46,2.605,4])
        cross_section_data_N2 = np.array([1.426,1.464,1.49,1.55,1.62,1.718,1.81,1.908,2.00,2.062,2.131,2.19,2.342,2.55,2.729,2.85,3.12,3.4,3.85,4.33,4.72,5.1,5.41,5.69,5.95,6.45,7.1,7.59,7.9,8.5,9.0,9.7,10.16,10.65,10.87,11.0,11.03,11.07,11.12,17.4,18.03,16.65,12.38,10.9])
        N2_cross_section = np.interp(Te,ev_data,cross_section_data_N2)
        return N2_cross_section*1e-20
    
    ######################################################################
    ######################################################################
    # Aerodynamics #######################################################
    ######################################################################
    ######################################################################  
    
    def braking_temperature(self,T_inf,Mach,gamma=1.4):
        # T_inf - scalar - free stream temperature (K)
        # Mach - scalar - free stream mach number
        # gamma - scalar - ratio of specific heats (set to 1.4 for air unless specified otherwise)
        return T_inf*(1+(gamma-1)/2 * Mach**2)
    
    def speed_of_sound(self,Temp,R_specific=287,gamma=1.4):
        # Temp - scalar - gas temperature (K)
        # R_specific - scalar - specific gas constant (automatically set to 287) (J/Kg*K)
        # gamma - scalar - ratio of specific heats (automatically set to 1.4)
        return np.sqrt(gamma*R_specific*Temp)
    
    def normal_shock_static_pressure_ratio(self,Mach,gamma=1.4):
        # Mach - scalar - free stream mach number
        # gamma - scalar - ratio of specific heats (automatically set to 1.4)
        return (2*gamma*Mach**2 - (gamma-1))/(gamma+1)
    
    def normal_shock_static_pressure_ratio(self,Mach,gamma=1.4):
        # Mach - scalar - free stream mach number
        # gamma - scalar - ratio of specific heats (automatically set to 1.4)
        return (((gamma+1)*Mach**2)/((gamma-1)*Mach**2+2))**(gamma/(gamma-1)) * ((gamma*1)/(2*gamma*Mach**2-(gamma-1)))**(1/(gamma-1))
    
    def normal_shock_density_ratio(self,Mach,gamma=1.4):
        # Mach - scalar - free stream mach number
        # gamma - scalar - ratio of specific heats (automatically set to 1.4)
        return ((gamma+1)*Mach**2)/((gamma-1)*Mach**2+2)
    
    def normal_shock_temperature_ratio(self,Mach,gamma=1.4):
        # Mach - scalar - free stream mach number
        # gamma - scalar - ratio of specific heats (automatically set to 1.4)
        return ((2*gamma*Mach**2-(gamma-1))*((gamma-1)*Mach**2+2))/(((gamma+1)**2 * Mach**2))
    
    
    ######################################################################
    ######################################################################
    # Mathematics ########################################################
    ######################################################################
    ######################################################################          

    def distance(self,Point1,Point2):
        # Point1 - vector/tuple - initial point (x,y,z)
        # Point2 - vector/tuple - final point (x,y,z)
        x1,y1,z1 = Point1
        x2,y2,z2 = Point2
        return np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

    def area_regular_polygon(self,sides,side_length):
        # sides - scalar - number of sides
        # side_length - scalar - length of single side
        R = side_length/(2*np.sin(np.pi/sides))
        return 0.5*sides*R**2 * np.sin(2*np.pi/sides)
    
    def arc_length(self,radi,angle):
        # radi - scalar - radius os the circle
        # angle - scalar - angle of arc (radians)
        return 2*radi*np.sin(angle/2)
    
    def vector_angle(self,r1,r2):
        # r1 - vector/tuple - vector line (x,y,z)
        # r2 - vector/tuple - vector line (x,y,z) 
        x1,y1,z1 = r1
        x2,y2,z2 = r2
        return np.arccos((x1*x2+y1*y2+z1*z2)/(np.sqrt(x1**2 + y1**2 + z1**2)*np.sqrt(x2**2 + y2**2 + z2**2)))
    
    def scalar_product(self,r1,r2):
        # r1 - vector/tuple - vector line (x,y,z)
        # r2 - vector/tuple - vector line (x,y,z) 
        x1,y1,z1 = r1
        x2,y2,z2 = r2
        return x1*x2+y1*y2+z1*z2
    
    def vector_product(self,r1,r2):
        # r1 - vector/tuple - vector line (x,y,z)
        # r2 - vector/tuple - vector line (x,y,z) 
        x1,y1,z1 = r1
        x2,y2,z2 = r2
        x_component = y1*z2-z1*y2
        y_component = -1*(x1*z2-z1*x2)
        z_component = x1*y2-y1*x2
        return np.array([x_component,y_component,z_component])
    
    def triangle_centroid(self,A_point,B_point,C_point):
        # A_point - vector/tuple - triangle vertex A (x,y)
        # B_point - vector/tuple - triangle vertex B (x,y)
        # C_point - vector/tuple - triangle vertex C (x,y)
        Ax,Ay = A_point
        Bx,By = B_point
        Cx,Cy = C_point
        X0 = (Ax+Bx+Cx)/3
        Y0 = (Ay+By+Cy)/3
        return np.array([X0,Y0])
    
    def magnitude(self,r1):
        # r1 - vector/tuple - vector line (x,y,z)
        x1,y1,z1 = r1
        return np.sqrt(x1**2+y1**2+z1**2)
    
    def __ccw(self,A,B,C):
        '''
        A,B,C - tuple - coordinate points (x,y)
        '''
        ax,ay=A
        bx,by=B
        cx,cy=C
        return (cy-ay)*(bx-ax) > (by-ay)*(cx-ax)
    
    def intersection(self,seg_1,seg_2):
        '''
        seg_1, seg_2 - list/tuple - line segments [(x1,y1),(x2,y2)]
        '''
        s1p1 = seg_1[0]
        s1p2 = seg_1[1]
        s2p1 = seg_2[0]
        s2p2 = seg_2[1]
        return self.__ccw(s1p1,s2p1,s2p2) != self.__ccw(s1p2,s2p1,s2p2) and self.__ccw(s1p1,s1p2,s2p1) != self.__ccw(s1p1,s1p2,s2p2)


    