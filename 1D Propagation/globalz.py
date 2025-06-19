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
        Calculate the propagation constant (β) in a medium with a given refractive index.
        Parameters:
            f_0 (float): Electromagnetic frequency in Hz
            refra (float): Refractive index of the medium
        Returns:
            float: Propagation constant (rad/m)
        '''
        W_0 = f_0*2*np.pi
        return W_0/self.c * refra
    
    def critical_angle(self,n1,n2):
        '''
        Compute the critical angle for total internal reflection.
        Parameters:
            n1 (float): Refractive index of incident medium
            n2 (float): Refractive index of transmitting medium
        Returns:
            float: Critical angle in radians
        '''
        return np.arcsin(n2/n1)
    
    def beer_lambert(self,f_0,refra,del_x):
        '''
        Compute attenuation of light using Beer-Lambert law.
        Parameters:
            f_0 (float): Frequency of light (Hz)
            refra (complex): Complex refractive index
            del_x (float): Distance light travels in the medium (m)
        Returns:
            float: Attenuated intensity ratio
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
        '''
        Compute the plasma frequency.
        Parameters:
            Ne (float): Electron number density (m⁻³)
        Returns:
            float: Plasma frequency in rad/s
        '''
        return np.sqrt((Ne*self.q**2)/(self.m_e*self.esp_0))
    
    def critical_electron_density(self,f_0):
        '''
        Compute the critical electron density for wave propagation.
        Parameters:
            f_0 (float): Frequency of EM wave (Hz)
        Returns:
            float: Critical electron density (m⁻³)
        '''
        W_0 = f_0*2*np.pi
        return W_0**2 * self.m_e*self.esp_0/(self.q**2)

    def gyrofrequency(self,B):
        '''
        Calculate the electron gyrofrequency.
        Parameters:
            B (float): Magnetic field strength (T)
        Returns:
            float: Gyrofrequency (rad/s)
        '''
        return self.q*B/(self.m_e)
    
    def electron_thermal_velocity(self,Te):
        '''
        Compute the most probable thermal velocity of electrons.
        Parameters:
            Te (float): Electron temperature (K)
        Returns:
            float: Thermal velocity (m/s)
        '''
        return np.sqrt(8*self.Kb*Te/(self.m_e*np.pi))
    
    def debye_length(self,Te,Ne):
        '''
        Compute Debye length in a plasma.
        Parameters:
            Te (float): Electron temperature (K)
            Ne (float): Electron number density (m⁻³)
        Returns:
            float: Debye length (m)
        '''
        return np.sqrt(self.esp_0*self.Kb*Te/(Ne*self.q**2))
    
    def ExB_drift(self,E_field,B_field):
        '''
        Compute E × B drift velocity vector.
        Parameters:
            E_field (ndarray): Electric field vector (V/m)
            B_field (ndarray): Magnetic field vector (T)
        Returns:
            ndarray: Drift velocity vector (m/s)
        '''
        return np.cross(E_field,B_field)/np.abs(B_field)
    
    def Hall_Parameter(self,B_mag,nu):
        '''
        Compute the Hall parameter.
        Parameters:
            B_mag (float): Magnetic field strength (T)
            nu (float): Collision frequency (s⁻¹)
        Returns:
            float: Hall parameter (dimensionless)
        '''
        return self.q*B_mag/(self.m_e*nu)
    
    
    
    ######################################################################
    ######################################################################
    # PLASMA KINETICS ####################################################
    ######################################################################
    ######################################################################

    def coulomb_collision_frequency(self,Ne,Te):
        """
        Calculate the electron-electron Coulomb collision frequency (ν_ee).

        Parameters:
        Ne : float
            Electron number density in m^-3.
        Te : float
            Electron temperature in Kelvin.

        Returns:
        float
            Collision frequency in s^-1.

        Formula source:
        The formula includes a Coulomb logarithm and is based on plasma kinetic theory.
        The Coulomb logarithm accounts for the range of impact parameters in electron-electron collisions.

        Notes:
        The logarithmic term (coul_log) involves fundamental constants and approximates the average number of
        small-angle deflections in a Coulomb interaction between charged particles.
        """
        coul_log = 0.5*np.log(1.0+0.25*self.esp_0*self.Kb*Te/(np.pi*self.q**2*Ne*(self.q**4/(16*np.pi**2*self.esp_0**2*self.Kb**2*Te**2)+self.h**2/(2*np.pi*self.Kb*self.m_e*Te))))
        return 3.633e-6*coul_log*Ne*Te**(-1.5)
    
    def electron_mobility(self,nu):
        """
        Calculate electron mobility (μ_e) from the collision frequency.

        Parameters:
        nu : float
            Electron collision frequency in s^-1.

        Returns:
        float
            Electron mobility in m^2/(V·s).

        Notes:
        Electron mobility is defined as the average drift velocity per unit electric field.
        This formulation is derived from classical transport theory:
        μ_e = q / (m_e * ν)
        """
        return self.q/(self.m_e*nu)
    
    def electron_diffusion_coefficient(self,Temp,nu):
        """
        Calculate the electron diffusion coefficient (D_e).

        Parameters:
        Temp : float
            Gas temperature in Kelvin.
        nu : float
            Electron collision frequency in s^-1.

        Returns:
        float
            Electron diffusion coefficient in m^2/s.

        Notes:
        Derived from the Einstein relation for diffusion:
        D_e = k_B * T / (m_e * ν)
        It assumes classical kinetic theory where random thermal motion drives diffusion.
        """
        return self.Kb*Temp/(self.m_e*nu)
    
    def bohm_diffusion_coefficient(self,B_mag,Te):
        """
        Calculate the Bohm diffusion coefficient (D_B).

        Parameters:
        B_mag : float
            Magnetic field strength in Tesla.
        Te : float
            Electron temperature in Kelvin.

        Returns:
        float
            Bohm diffusion coefficient in m^2/s.

        Notes:
        Bohm diffusion is an empirical model used for cross-field (perpendicular to magnetic field)
        transport in magnetized plasmas, typically much larger than classical diffusion.

        Formula:
        D_B = k_B * Te / (16 * q * B)
        """
        return self.Kb*Te/(16*self.q*B_mag)

    ######################################################################
    ######################################################################
    # PLASMA OPTICS ######################################################
    ######################################################################
    ######################################################################

    def AH_parallel(self,Ne,B_mag,f_0,nu):
        """
        Calculate the squared refractive index (n^2) for electromagnetic wave propagation 
        perpendicular to the magnetic field using the Appleton-Hartree formula (simplified case).

        Parameters:
        Ne : float
            Electron number density (m^-3).
        B_mag : float
            Magnetic field magnitude (T).
        f_0 : float
            Electromagnetic wave frequency (Hz).
        nu : float
            Effective electron-neutral collision frequency (s^-1).

        Returns:
        complex
            Squared refractive index n^2 (complex-valued to include absorption).
        
        Notes:
        For propagation perpendicular to the magnetic field, the Appleton-Hartree expression simplifies,
        but still involves cyclotron and plasma resonances.
        """
        omega_p = self.plasma_frequency(Ne)
        omega_c = self.gyrofrequency(B_mag)
        frequency = f_0*2*np.pi
        n2 = 1-(omega_p/frequency)**2/(1-1j*(nu/frequency)-(omega_c/frequency)**2/(1-1j*(nu/frequency)-(omega_p/frequency)**2))
        return n2

    def AH_perpendicular(self,Ne,B_mag,f_0,nu):
        """
        Calculate the squared refractive index (n^2) for electromagnetic wave propagation 
        parallel to the magnetic field using a simplified Appleton-Hartree relation.

        Parameters:
        Ne : float
            Electron number density (m^-3).
        B_mag : float
            Magnetic field magnitude (T).
        f_0 : float
            Electromagnetic wave frequency (Hz).
        nu : float
            Electron-neutral collision frequency (s^-1).

        Returns:
        complex
            Squared refractive index n^2 (complex-valued).
        
        Notes:
        In the parallel case, the interaction with the magnetic field is less complex,
        leading to a simplified dispersion relation.
        """
        omega_p = self.plasma_frequency(Ne)
        omega_c = self.gyrofrequency(B_mag)
        frequency = f_0*2*np.pi
        n2 = 1-(omega_p/frequency)**2/(1-1j*(nu/frequency)-(omega_c/frequency))
        return n2

    def Appleton_Hartree(self,Ne,B_mag,f_0,nu,angle,polar = 'rhc'):
        """
        General Appleton-Hartree formula for arbitrary propagation angle and circular polarization.

        Parameters:
        Ne : float
            Electron number density (m^-3).
        B_mag : float
            Magnetic field magnitude (T).
        f_0 : float
            Wave frequency (Hz).
        nu : float
            Electron-neutral collision frequency (s^-1).
        angle : float
            Propagation angle with respect to magnetic field (radians).
        polar : str
            'rhc' for right-hand circular, 'lhc' for left-hand circular polarization.

        Returns:
        complex
            Squared refractive index n^2 (complex-valued).

        Notes:
        This is the full Appleton-Hartree formula. It captures the complex interplay between
        the plasma oscillation, cyclotron motion, and collisions for any angle of incidence.
        """
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
        """
        Check for cold plasma cutoff condition where wave propagation is not supported.

        Parameters:
        Ne : float
            Electron number density (m^-3).
        B_mag : float
            Magnetic field magnitude (T).
        f_0 : float
            Wave frequency (Hz).

        Returns:
        bool
            True if cutoff condition is met, i.e., wave cannot propagate; False otherwise.

        Notes:
        Cold plasma cutoff occurs when the plasma frequency equals the upper or lower hybrid resonance.
        """
        W_p = self.plasma_frequency(Ne)
        W_c = self.gyrofrequency(B_mag)
        W_0 = f_0*2*np.pi
        if W_p**2-(W_c+W_0)<=0.00000001 or W_p**2-(W_0-W_c)<=0.00000001: 
            return True
        else:
            return False
        
    def cold_resonance(self,Ne,B_mag,f_0):
        """
        Check for cold plasma resonance condition.

        Parameters:
        Ne : float
            Electron number density (m^-3).
        B_mag : float
            Magnetic field magnitude (T).
        f_0 : float
            Wave frequency (Hz).

        Returns:
        bool
            True if resonance condition is met, i.e., enhanced wave-particle interaction occurs.

        Notes:
        Resonance occurs when the driving frequency matches a natural frequency of the plasma system.
        """
        W_p = self.plasma_frequency(Ne)
        W_c = self.gyrofrequency(B_mag)
        W_0 = f_0*2*np.pi
        return W_p**2-(W_0**2-W_c**2)<=0.00000001

    def air_to_plasma_transmission(self,f_0,n_plasma):
        """
        Calculate transmission coefficient at the air-plasma interface.

        Parameters:
        f_0 : float
            Frequency of the electromagnetic wave (Hz).
        n_plasma : float
            Refractive index of the plasma (can be complex).

        Returns:
        complex
            Transmission coefficient (can be complex, indicating phase shift and absorption).

        Notes:
        Uses Fresnel relations assuming normal incidence.
        """
        K_0 = self.propagation_constant(f_0,1.0)
        K_1 = self.propagation_constant(f_0,n_plasma)
        return 2*K_0/(K_0+K_1)
    
    def air_to_plasma_reflection(self,f_0,n_plasma):
        """
        Calculate reflection coefficient at the air-plasma interface.

        Parameters:
        f_0 : float
            Frequency of the electromagnetic wave (Hz).
        n_plasma : float
            Refractive index of the plasma.

        Returns:
        complex
            Reflection coefficient (complex for phase/amplitude shift).

        Notes:
        Uses Fresnel equations assuming normal incidence.
        """
        K_0 = self.propagation_constant(f_0,1.0)
        K_1 = self.propagation_constant(f_0,n_plasma)
        return (K_0-K_1)/(K_0+K_1)

    ######################################################################
    ######################################################################
    # KINETICS ###########################################################
    ######################################################################
    ######################################################################

    def N2_cross_section_low_ev(self,Te):
        """
        Estimate the electron impact cross-section of molecular nitrogen (N₂) 
        at low electron temperatures (< 4 eV).

        Parameters:
        Te : float
            Electron temperature in Kelvin. This value is converted to eV using Te / 11600.

        Returns:
        float
            Estimated electron impact cross-section for N₂ (in m²).
        
        Notes:
        - Uses linear interpolation of experimental or tabulated cross-section data.
        - Cross-section values are typically given in units of 10⁻²⁰ m².
        - Relevant for modeling electron impact excitation/ionization in low-temperature plasmas.
        """
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
        """
        Calculate the stagnation (braking) temperature from free stream conditions.

        Parameters:
        T_inf : float
            Free stream static temperature (K).
        Mach : float
            Free stream Mach number.
        gamma : float
            Ratio of specific heats (default = 1.4 for air).

        Returns:
        float
            Braking (stagnation) temperature (K).
        
        Notes:
        This is derived from the isentropic flow relations for compressible flow.
        """
        return T_inf*(1+(gamma-1)/2 * Mach**2)
    
    def speed_of_sound(self,Temp,R_specific=287,gamma=1.4):
        """
        Compute the speed of sound in a gas.

        Parameters:
        Temp : float
            Temperature of the gas (K).
        R_specific : float
            Specific gas constant (J/kg·K), default is for air.
        gamma : float
            Ratio of specific heats (default = 1.4).

        Returns:
        float
            Speed of sound (m/s).
        
        Notes:
        Derived from a = sqrt(γ * R * T), valid for ideal gases.
        """
        return np.sqrt(gamma*R_specific*Temp)
    
    def normal_shock_static_pressure_ratio(self,Mach,gamma=1.4):
        """
        Calculate the static pressure ratio (p2/p1) across a normal shock.

        Parameters:
        Mach : float
            Upstream Mach number.
        gamma : float
            Ratio of specific heats (default = 1.4).

        Returns:
        float
            Static pressure ratio across the shock.

        Notes:
        This is derived from normal shock relations in compressible flow.
        """
        return (2*gamma*Mach**2 - (gamma-1))/(gamma+1)
    
    def normal_shock_static_pressure_ratio(self,Mach,gamma=1.4):
        """
        Calculate the total pressure ratio (p02/p01) across a normal shock.

        Parameters:
        Mach : float
            Upstream Mach number.
        gamma : float
            Ratio of specific heats (default = 1.4).

        Returns:
        float
            Total pressure ratio across the shock.

        Notes:
        - This reflects the irreversible loss in total pressure due to the shock.
        - Important for estimating performance loss in compressible systems.
        """
        return (((gamma+1)*Mach**2)/((gamma-1)*Mach**2+2))**(gamma/(gamma-1)) * ((gamma*1)/(2*gamma*Mach**2-(gamma-1)))**(1/(gamma-1))
    
    def normal_shock_density_ratio(self,Mach,gamma=1.4):
        """
        Calculate the density ratio (ρ2/ρ1) across a normal shock.

        Parameters:
        Mach : float
            Upstream Mach number.
        gamma : float
            Ratio of specific heats (default = 1.4).

        Returns:
        float
            Density ratio across the shock.
        
        Notes:
        Useful for continuity and mass flow analysis in compressible systems.
        """
        return ((gamma+1)*Mach**2)/((gamma-1)*Mach**2+2)
    
    def normal_shock_temperature_ratio(self,Mach,gamma=1.4):
        """
        Calculate the temperature ratio (T2/T1) across a normal shock.

        Parameters:
        Mach : float
            Upstream Mach number.
        gamma : float
            Ratio of specific heats (default = 1.4).

        Returns:
        float
            Temperature ratio across the shock.

        Notes:
        Helps in estimating thermal loads and post-shock conditions.
        """
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


    