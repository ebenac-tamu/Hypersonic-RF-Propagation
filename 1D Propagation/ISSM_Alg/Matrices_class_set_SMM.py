import numpy as np

class scattering_matrix_set:
    def __init__(self):
        """
        Constructor method for initializing class attributes.

        Attributes:
        ----------
        self.n : int
            Initialized to -2. 

        self.p : int
            Initialized to -1. 
        """
        self.n = -2
        self.p = -1
    
    def __incident_S_matrix(self,k):
        """
        Constructs the incident scattering matrix (S-matrix) for a two-layer interface
        in a 1D stratified medium.

        Parameters
        ----------
        k : list or array-like of complex
            A list or array of complex wave numbers (propagation constants) where:
            - k[0] is the wave number in the incident (first) medium
            - k[1] is the wave number in the first layer of the stratified medium

        Returns
        -------
        numpy.ndarray (2x2 complex)
            The scaled scattering matrix S_0 for the interface between the incident medium
            and the first layer. This matrix maps forward and backward wave amplitudes
            across the boundary.
        """
        const = 1/(2*k[1])
        S_0_array = np.array([[k[1]+k[0], k[1]-k[0]],
                            [k[1]-k[0], k[1]+k[0]]],dtype=complex)
        return const*S_0_array 
    
    def __Scattering_Matrix_n(self,z,k,d):
        """
        Computes the local scattering matrix between layer (z-1) and layer z in a stratified medium,
        based on the continuity of electric and magnetic fields at the interface.

        Parameters
        ----------
        z : int
            Index of the current layer (starting from 1).
        k : list or np.ndarray of complex
            Complex wave numbers (propagation constants) for each layer.
        d : list or np.ndarray of float
            Cumulative distances (thickness positions) at each layer boundary.

        Returns
        -------
        smn : 2x2 complex np.ndarray
            Local scattering matrix between layer z-1 and z, which connects the wave amplitudes
            across the interface using field-matching boundary conditions.
        """
        matrix_1 = np.array([[np.exp(-1j*k[z]*(d[z-1])), np.exp(1j*k[z]*(d[z-1]))],
                             [k[z]*np.exp(-1j*k[z]*(d[z-1])), -k[z]*np.exp(1j*k[z]*(d[z-1]))]], dtype=complex)
        #print(matrix_1)
        matrix_2 = np.array([[np.exp(-1j*k[z-1]*(d[z-1])), np.exp(1j*k[z-1]*(d[z-1]))],
                             [k[z-1]*np.exp(-1j*k[z-1]*(d[z-1])), -1*k[z-1]*np.exp(1j*k[z-1]*(d[z-1]))]],dtype=complex)
        matrix_1_inv = np.linalg.inv(matrix_1)
        #print(k[n-1]*(np.exp(1j*k[n-1]*(d[n]-d[n-1]))))
        smn = np.matmul(matrix_1_inv,matrix_2)
        #print(smn)
        return smn
    
    def __boundary_ISMM(self,k,d):
        """
        Computes the boundary condition vector (Vp) for the Incident Side Scattering Matrix Method (ISMM).
        
        This vector represents the wave field at the left boundary (incident interface) of the layered medium,
        where a unit-amplitude wave is assumed to be incident from the left.

        Parameters
        ----------
        k : list or np.ndarray of complex
            Complex propagation constants (wave numbers) for each layer.
        d : list or np.ndarray of float
            Cumulative spatial positions of layer boundaries (thickness values).

        Returns
        -------
        Vp : 2x1 complex np.ndarray
            Boundary condition vector that initializes the recursive S-matrix algorithm.
        """
        #Sp = self.__Scattering_Matrix_n(-1,k,d)
        #Sp_inv = np.linalg.inv(Sp)
        #Vp = np.matmul(Sp_inv,np.array([[1],[0]]))
        constant = 1/(2*k[self.n])
        array = np.array([[(k[self.n]+k[self.p])*np.exp(1j*(k[self.n]-k[self.p])*(d[self.p]))],
                          [(k[self.n]-k[self.p])*np.exp(-1j*(k[self.n]+k[self.p])*(d[self.p]))]])
        Vp = constant*array
        return Vp
    
    def __global_scattering_matrix(self,k,d):
        """
        Constructs the global scattering matrix S_g_matrix that relates the total wave fields
        from the incident side to the boundary of the last layer.

        This function implements the core of the Scattering Matrix Method (SMM), which recursively
        computes how waves reflect and transmit through a multilayer structure.

        Parameters
        ----------
        k : list or np.ndarray of complex
            Complex wave numbers (propagation constants) for each layer.
        d : list or np.ndarray of float
            Cumulative positions of the interfaces between layers. The length of d determines the number of layers.

        Returns
        -------
        S_g_matrix : 2x2 complex np.ndarray
            Global scattering matrix that represents total transmission and reflection effects
            of all layers in the structure from the incident interface to the end.
        """
        '''
        S1 = self.__incident_S_matrix(k)
        Sn = self.__Scattering_Matrix_n(len(d)-2,k,d)
        print(S1)
        print(self.__Scattering_Matrix_n(0,k,d))
        for i in range(len(d)-3,2,-1):
            #print(i)
            Sn = np.matmul(Sn,self.__Scattering_Matrix_n(i,k,d))
        S_g_matrix = np.matmul(Sn,S1)
        '''
        incidence = self.__incident_S_matrix(k)
        #print(incidence)
        S_g_matrix = np.matmul(self.__Scattering_Matrix_n(2,k,d),incidence)
        for i in range(len(d)-3):
            S_g_matrix = np.matmul(self.__Scattering_Matrix_n(i+3,k,d),S_g_matrix)
        #print(S_g_matrix)
        
        return S_g_matrix
    
    def Total_coefficients(self,k,d):
        """
        Computes the total complex wave coefficients A and D using the global scattering matrix and
        boundary conditions for a layered structure.

        These coefficients represent the amplitudes of the forward (A) and backward (D) propagating
        waves in the first and last layers, respectively.

        Parameters
        ----------
        k : list of complex
            Complex propagation constants for each layer.
        d : list of float
            Layer boundary positions; must be of length equal to number of layers + 1.

        Returns
        -------
        coeffs : np.ndarray (complex)
            Array of shape (2,) containing [A, D]:
            - A: Incident wave coefficient (input layer)
            - D: Transmitted/reflected wave coefficient (output layer)
        """
        S_g = self.__global_scattering_matrix(k,d)
        #print(S_g)
        S_g1_col = S_g[:,0]
        #print(S_g1_col)
        S_g2_col = S_g[:,1]
        #print(S_g2_col)
        V_p_col = self.__boundary_ISMM(k,d)
        #print(V_p_col)
        #ad = S_g2_col-V_p_col
        #ad = -1*np.linalg.inv(ad)@S_g1_col
        #print(np.abs(ad)**2)
        #ad1 = (S_g2_col[0]-V_p_col[0])
        #ad2 = (S_g2_col[1]-V_p_col[1])
        #print(ad1)
        #print(ad2)
        #A = -1*(1/ad1)*S_g1_col[0]
        #D = -1*(1/ad2)*S_g1_col[1]

        M = np.array([[S_g2_col[0],-1*V_p_col[0][0]],
                      [S_g2_col[1],-1*V_p_col[1][0]]])
        #print(np.linalg.inv(M))
        rt = np.dot(-1*np.linalg.inv(M),S_g1_col)
        #print(rt)
        #print(np.abs(rt)**2)
        #print(np.abs(ad1))
        #print(np.abs(ad2))
        A = (S_g1_col[0]*V_p_col[1][0]-S_g1_col[1]*V_p_col[0][0])/(S_g2_col[1]*V_p_col[0][0]-S_g2_col[0]*V_p_col[1][0])
        #print(A)
        #print(np.abs(A))
        D = (S_g1_col[0]+S_g2_col[0]*A)/V_p_col[0][0]
        #print(D)
        #print(np.abs(D))
        #sudo_inv = np.linalg.pinv(S_g1_col-V_p_col , rcond=1e-10)
        #coeffs = np.array([[A],[D]],dtype=complex)
        #print('---')
        #print(A)
        #print(D)
        #print(rt)
        #print('---')
        return np.array([A,D])
    
    def partial_coefficients(self,k,d,A,D):
        """
        Computes the local wave coefficients (B and C) at each interface within a layered medium
        using the forward-propagating (A) and backward-propagating (D) total wave coefficients.

        These partial coefficients represent the right-going (B) and left-going (C) wave amplitudes
        at each layer interface inside the structure.

        Parameters
        ----------
        k : list of complex
            Propagation constants for each layer.
        d : list of float
            Layer interface locations (length = number of layers + 1).
        A : complex
            Forward-propagating coefficient from total solution (input wave).
        D : complex
            Backward-propagating coefficient from total solution (output/reflected wave).

        Returns
        -------
        B_mag : np.ndarray of complex
            Right-going wave amplitudes at each layer interface.
        C_mag : np.ndarray of complex
            Left-going wave amplitudes at each layer interface.
        B_phase : np.ndarray of complex
            Phase of B waves at each layer interface.
        C_phase : np.ndarray of complex
            Phase of C waves at each layer interface.
        """
        B_mag = np.zeros((len(d)),dtype=complex)
        C_mag = np.zeros((len(d)),dtype=complex)
        B_phase = np.zeros((len(d)),dtype=complex)
        C_phase = np.zeros((len(d)),dtype=complex)
        incident_vector = np.array([[A],[1]],dtype=complex)
        S_1_matrix = self.__incident_S_matrix(k)
        #print(incident_vector)
        #print(S_1_matrix)
        BC_vector = S_1_matrix@incident_vector
        B_mag[0] = BC_vector[0][0]
        C_mag[0] = BC_vector[1][0]
        B_phase[0] = np.arctan(np.imag(BC_vector[0][0])/np.real(BC_vector[0][0]))
        C_phase[0] = np.arctan(np.imag(BC_vector[1][0])/np.real(BC_vector[1][0]))
        for m in range(1,len(B_mag)):
            S_m = self.__Scattering_Matrix_n(m,k,d)
            BC_vector = np.matmul(S_m,BC_vector)

            B_mag[m] = BC_vector[0][0]
            C_mag[m] = BC_vector[1][0]
            B_phase[m] = np.arctan(np.imag(BC_vector[0][0])/np.real(BC_vector[0][0]))
            C_phase[m] = np.arctan(np.imag(BC_vector[1][0])/np.real(BC_vector[1][0]))
        V_p = self.__boundary_ISMM(k,d)
        BC_vector_last = V_p*D#np.array([[D],[0]])
        B_mag[-1] = BC_vector_last[0][0]
        C_mag[-1] = BC_vector_last[1][0]
        return B_mag,C_mag,B_phase,C_phase