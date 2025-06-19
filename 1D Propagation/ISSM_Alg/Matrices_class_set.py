import numpy as np

##########################################################
#################### Not in Use ##########################
##########################################################

class scattering_matrix_set:
    def __init__(self):
        self.n = -2
        self.p = -1
    
    def __incident_S_matrix(self,k):
        const = 1/(2*k[1])
        S_0_array = np.array([[k[1]+k[0], k[1]-k[0]],
                            [k[1]-k[0], k[1]+k[0]]],dtype=complex)
        return const*S_0_array 
    
    def __Scattering_Matrix_n(self,n,k,d):
        matrix_1 = np.array([[1.0, 1.0],
                             [k[n], -k[n]]], dtype=complex)
        #print(np.linalg.cond(matrix_1))
        matrix_2 = np.array([[np.exp(-1j*k[n-1]*(d[n-1]-d[n-2])), np.exp(1j*k[n-1]*(d[n-1]-d[n-2]))],
                             [k[n-1]*np.exp(-1j*k[n-1]*(d[n-1]-d[n-2])), -1*k[n-1]*np.exp(1j*k[n-1]*(d[n-1]-d[n-2]))]],dtype=complex)
        matrix_1_inv = np.linalg.inv(matrix_1)
        #print(k[n-1]*(np.exp(1j*k[n-1]*(d[n]-d[n-1]))))
        smn = np.matmul(matrix_1_inv,matrix_2)
        return smn
    
    def __boundary_ISMM(self,k,d):
        #Sp = self.__Scattering_Matrix_n(-1,k,d)
        #Sp_inv = np.linalg.inv(Sp)
        #Vp = np.matmul(Sp_inv,np.array([[1],[0]]))
        constant = 1/(2*k[self.n])
        array = np.array([[(k[self.n]+k[self.p])*np.exp(1j*k[self.n]*(d[self.p]-d[self.n]))],
                          [(k[self.n]-k[self.p])*np.exp(-1j*k[self.n]*(d[self.p]-d[self.n]))]])
        Vp = constant*array
        return Vp
    
    def __global_scattering_matrix(self,k,d):
        '''
        S1 = self.__incident_S_matrix(k)
        Sn = self.__Scattering_Matrix_n(len(d)-2,k,d)
        for i in range(len(d)-3,2,-1):
            #print(i)
            Sn = np.matmul(Sn,self.__Scattering_Matrix_n(i,k,d))
        S_g_matrix = np.matmul(Sn,S1)
        
        incidence = self.__incident_S_matrix(k)
        S_g_matrix = np.matmul(self.__Scattering_Matrix_n(1,k,d),incidence)
        for i in range(len(d)-3):
            S_g_matrix = np.matmul(self.__Scattering_Matrix_n(i+1,k,d),S_g_matrix)
        #print(S_g_matrix)
        '''
        incidence = self.__incident_S_matrix(k)
        #print(incidence)
        S_g_matrix = np.matmul(self.__Scattering_Matrix_n(2,k,d),incidence)
        for i in range(len(d)-3):
            S_g_matrix = np.matmul(self.__Scattering_Matrix_n(i+3,k,d),S_g_matrix)
        #print(S_g_matrix)
        return S_g_matrix
    
    def Total_coefficients(self,k,d):
        S_g = self.__global_scattering_matrix(k,d)
        #print(np.linalg.cond(S_g))
        S_g1_col = S_g[:,0]
        #print(S_g1_col[0])
        S_g2_col = S_g[:,1]
        #print(S_g2_col[0])
        V_p_col = self.__boundary_ISMM(k,d)
        #print(V_p_col[0][0])
        #ad = S_g2_col-V_p_col
        #ad = -1*np.linalg.inv(ad)@S_g1_col
        #print(np.abs(ad)**2)
        ad1 = (S_g2_col[0]-V_p_col[0])
        ad2 = (S_g2_col[1]-V_p_col[1])
        #print(ad1)
        #print(ad2)
        A = -1*(1/ad1)*S_g1_col[0]
        D = -1*(1/ad2)*S_g1_col[1]

        M = np.array([[S_g2_col[0],-1*V_p_col[0][0]],
                      [S_g2_col[1],-1*V_p_col[1][0]]])
        rt = -1*np.linalg.inv(M)@S_g1_col.T
        #print(rt)
        #print(np.abs(rt)**2)
        #print(np.abs(ad1))
        #print(np.abs(ad2))
        #A = (S_g1_col[0]*V_p_col[1]-S_g1_col[1]*V_p_col[0])/(S_g2_col[1]*V_p_col[0]-S_g2_col[0]*V_p_col[1])
        #print(A)
        #print(np.abs(A))
        #D = (S_g1_col[0]+S_g2_col[0]*A)/V_p_col[0]
        #print(D)
        #print(np.abs(D))
        #sudo_inv = np.linalg.pinv(S_g1_col-V_p_col , rcond=1e-10)
        #coeffs = np.array([[A],[D]],dtype=complex)
        return rt
    
    def partial_coefficients(self,k,d,A,D):
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
        for m in range(1,len(B_mag)-1):
            S_m = self.__Scattering_Matrix_n(m,k,d)
            BC_vector = S_m@BC_vector

            B_mag[m] = BC_vector[0][0]
            C_mag[m] = BC_vector[1][0]
            B_phase[m] = np.arctan(np.imag(BC_vector[0][0])/np.real(BC_vector[0][0]))
            C_phase[m] = np.arctan(np.imag(BC_vector[1][0])/np.real(BC_vector[1][0]))
        V_p = self.__boundary_ISMM(k,d)
        BC_vector_last = V_p*D#np.array([[D],[0]])
        B_mag[-1] = BC_vector_last[0][0]
        C_mag[-1] = BC_vector_last[1][0]
        return B_mag,C_mag,B_phase,C_phase
