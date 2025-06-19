import numpy as np
import matplotlib.pyplot as plt

class RK4Solver:
    def __init__(self, dt,X_limit,Y_limit,wall,n_eq):
        """
        system: Function defining the system of equations. It should take (t, y) as inputs,
                where y is an array of 4 variables.
        dt: Time step for the integration.
        """
        self.dt = dt
        self.n_eq = n_eq
        #self.rf = rf
        self.xlim = X_limit
        self.ylim = Y_limit
        self.wall = wall
    
    # Eikonal
    def system(self,t, z):
        """Defines a system of 4 first-order ODEs."""
        p1,p2,x,y = z
        n_point,dndx,dndy = self.n_eq(x,y)
        dp1dt = np.real(dndx)/np.real(n_point)
        dp2dt = np.real(dndy)/np.real(n_point)
        #dp1dt = np.nan_to_num(dp1dt)
        #dp2dt = np.nan_to_num(dp2dt)
        dydt = [
            p1/np.real(n_point)**2,         # Example: dx/dt = p1/np.real(n_point)**2
            p2/np.real(n_point)**2,        # dy/dt = p2/np.real(n_point)**2
            dp1dt,         # dp1/dt = np.real(dndx)/np.real(n_point)
            dp2dt         # dp2/dt = np.real(dndy)/np.real(n_point)
        ]
        return dydt
    
    def step(self, t, y):
        """Performs a single RK4 step."""
        dt = self.dt
        k1 = np.array(self.system(t, y))
        k2 = np.array(self.system(t + dt/2, y + dt*k1/2))
        k3 = np.array(self.system(t + dt/2, y + dt*k2/2))
        k4 = np.array(self.system(t + dt, y + dt*k3))
        
        y_next = y + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
        return y_next
    
    def solve(self, y0, t0, t_end):
        """
        Solves the system from t0 to t_end with initial conditions y0.
        
        y0: Initial state (array of 4 elements)
        t0: Initial time
        t_end: Final time
        """
        t_values = [t0]
        y_values = [y0]
        
        t = t0
        y = np.array(y0)
        
        while t < t_end and y[2]>self.xlim[0] and y[2]<self.xlim[1] and y[3]>self.ylim[0] and y[3]<self.ylim[1] and y[3]>self.wall(y[2]):
            y = self.step(t, y)
            t += self.dt
            t_values.append(t)
            y_values.append(y)
        
        return np.array(t_values), np.array(y_values)

    def solve_eikonal(self,Initial_X,Initial_Y,Initial_angle):
        n_0,nul1,nul2 = self.n_eq(Initial_X,Initial_Y)
        y0 = np.array([n_0*np.cos(Initial_angle),n_0*np.sin(Initial_angle),Initial_X,Initial_Y])  # Initial conditions
        t0, t_end = 0, 100.0

        # Solve the system
        t_values, y_values = self.solve(y0, t0, t_end)
        return y_values



# Print some results
#print("Time values:", t_values[:5])
#print("Solution sample:", y_values[:5])

### Inputs for validation ###

dt = 0.01
n0 = 1000.0
X0 = 2
Y0 = 0
R = 1.0
X_limit = [0.0,5.0]
Y_limit = [-2,2]
wall = lambda x: -100

def n_luneberg(x,y):
    
    if ((x-X0)**2 +(y-Y0)**2)<=R:
        n = n0*np.sqrt(2-((x-X0)**2 + (y-Y0)**2)/(R**2))
        dndx = -1*(n0**2 * (x-X0))/(R**2 *n0*np.sqrt(2-((x-X0)**2 + (y-Y0)**2)/(R**2)))
        dndy = -1*(n0**2 * (y-Y0))/(R**2 *n0*np.sqrt(2-((x-X0)**2 + (y-Y0)**2)/(R**2)))
    else: 
        n = n0
        dndx = 0
        dndy = 0
    return n, dndx, dndy



def n_test(x,y):
    return 1.01,0.0,0.001

rta = RK4Solver(dt,X_limit,Y_limit,wall,n_luneberg)
output = rta.solve_eikonal(1.0,0,np.pi/3)

print(output[:,3])
plt.plot(output[:,2],output[:,3])
plt.show()