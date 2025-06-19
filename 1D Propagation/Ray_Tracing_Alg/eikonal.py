import numpy as np
from scipy.integrate import solve_ivp
#from Input_Data.Ray_tracing_data_management import import_data_refractive_index as rf
from globalz import global_variables as glbl
from typing import Protocol, Any

class Event(Protocol):
    terminal: bool = False
    direction: float = 0.0

    def __call__(self, t, y, *args):
        pass



class ray_tracing:
    def __init__(self,X_limit,Y_limit,rf):
        self.X_min,self.X_max = X_limit
        self.Y_min,self.Y_max = Y_limit
        self.rf = rf

    def eikonal_ray_tracing(self,Initial_X,Initial_Y,Initial_angle,step,data_rec=False):        
        
        initial_dir = np.array([Initial_X*np.cos(Initial_angle),Initial_Y*np.sin(Initial_angle),0])
        
        n,dndx,dndy = self.rf.refractive_index_set(Initial_X,Initial_Y,initial_dir)
        #print(n)
        if data_rec==True:
            #print('recording')
            n_real = [np.real(n)]
            n_imag = [np.imag(n)]
            X_data = [Initial_X]
            Y_data = [Initial_Y]
            #Ray_data = [(Initial_X,Initial_Y)]
        n = np.real(n)
        #print(n)
        #print()
        dndx = np.real(dndx)
        #print(dndx)
        #print()
        dndy = np.real(dndy)
        #print(dndy)
        #print()
        P1 = n*np.cos(Initial_angle)
        #print(P1)
        P2 = n*np.sin(Initial_angle)
        #print(P2)
        x1 = step*P1/(n**2)+Initial_X
        y1 = step*P2/ (n**2)+Initial_Y
        #print(X_min)
        #print(x1)
        #print(X_max)
        #print()
        #print(Y_min)
        #print(y1)
        #print(Y_max)
        r_diff = np.array([x1,y1,0]) - np.array([Initial_X,Initial_Y,0])
        while self.X_max>x1>self.X_min and self.Y_max>y1>=self.Y_min:
            n,dndx,dndy = self.rf.refractive_index_set(x1,y1,r_diff)
            if data_rec==True:
                n_real.append(np.real(n))
                n_imag.append(np.imag(n))
            n = np.real(n)
            dndx = np.real(dndx)
            dndy = np.real(dndy)
            
            P1 += step*dndx/n
            P2 += step*dndy/n
            
            r0 = np.array([x1,y1,0])
            x1 += step*P1/ (n**2)
            #print(X_min)
            #print(x1)
            #print(X_max)
            #print()
            #print(Y_min)
            #print(y1)
            #print(Y_max)
            y1 += step*P2/ (n**2)
            r_diff = np.array([x1,y1,0]) - r0
            if data_rec==True:
                #Ray_data.append(tuple((x1,y1)))
                #print(Ray_data)
                X_data.append(x1)
                Y_data.append(y1)
        #print(X_data)
        if data_rec==False:
            return x1,y1
        elif data_rec==True:
            #print(Ray_data)
            return X_data,Y_data, n_real, n_imag

    def __rk_45_eq(self,t,z):
        x,y,p1,p2 = z
        n_point,dndx,dndy = self.rf.refractive_index_set(x,y,[1,0])
        #print(t)
        #print(n_point)
        n_point = np.nan_to_num(n_point)
        #print(n_point**2)
        dxdt = p1/np.real(n_point**2)
        dydt = p2/np.real(n_point**2)
        dp1dt = np.real(dndx)/np.real(n_point)
        dp2dt = np.real(dndy)/np.real(n_point)
        dp1dt = np.nan_to_num(dp1dt)
        dp2dt = np.nan_to_num(dp2dt)
        #print(dp1dt)
        #print(dp2dt)
        return [dxdt,dydt,dp1dt,dp2dt]
    
    def EventDecorator(terminal: bool, direction: float):
        def decorator_event(func: Any) -> Event:
            func.terminal = terminal
            func.direction = direction
            return func

        return decorator_event
    
    @EventDecorator(terminal=True, direction=0.0)
    def __rk_45_y_event(self,t,vars):
        x,y,p1,p2 = vars
        if np.abs(y-self.Y_max)<0.0001 or y-self.Y_min<0:
            ret = 0
        else:
            ret = 1
        return ret

    @EventDecorator(terminal=True, direction=0.0)
    def __rk_45_x_event(self,t,vars):
        x,y,p1,p2 = vars
        if np.abs(x-self.X_min)<0.0001 or np.abs(x-self.X_max)<0.0001:
            ret = 0
        else:
            ret = 1
        return ret

    def __rk_45_p1_event(self,t,vars):
        x,y,p1,p2 = vars
        return 1
    
    def __rk_45_p2_event(self,t,vars):
        x,y,p1,p2 = vars
        return 1

    def eikonal_fast(self,Initial_X,Initial_Y,Initial_angle,data_rec=False):
        n, null1,nul2 = self.rf.refractive_index_set(Initial_X,Initial_Y,[1,0])
        initial_dir = np.array([np.real(n)*np.cos(Initial_angle),np.real(n)*np.sin(Initial_angle),0])
        t_span = [0,1000000.0]
        y0 = [Initial_X,Initial_Y,initial_dir[0],initial_dir[1]]
        #self.__rk_45_y_event.terminal = True
        event = [self.__rk_45_x_event,self.__rk_45_y_event,self.__rk_45_p1_event,self.__rk_45_p2_event]
        #print('pre-RK45')
        sol = solve_ivp(self.__rk_45_eq,t_span,y0,method='LSODA',max_step=1e-2,min_step=1e-4,events=event)
        #print('post-RK45')
        #print(sol.t[-1])
        #print(sol.y[0,-2])
        #print(sol.y[1,-2])
        #print(sol.y[0,-1])
        #print(sol.y[1,-1])
        x_val = sol.y[0,:]
        y_val = sol.y[1,:]
        if data_rec==False:
            return x_val[-1],y_val[-1]
        elif data_rec==True:
            return x_val,y_val
        
       