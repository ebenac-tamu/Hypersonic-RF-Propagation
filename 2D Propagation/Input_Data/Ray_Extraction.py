import pandas as pd

class Extract_Ray_Data:
    def __init__(self,rays,ray_x_data,ray_y_data):
        self.number_of_rays = rays
        self.x_rays = ray_x_data
        self.y_rays = ray_y_data
        
    def __ray_data_titles(self):
        self.titles = []
        for i in range(self.number_of_rays):
            n_ray = 'Ray Number '+str(i)
            self.titles.append(n_ray)
    
    def ray_data_extract(self):
        self.__ray_data_titles()
        total_ray_list = []
        for i in range(len(self.theta_range)):
            column = []
            #print(column)
            for j in range(len(self.x_rays[i])):
                column.append(tuple((self.x_rays[i][j],self.y_rays[i][j])))
            total_ray_list.append(column)
        #total_ray_list = np.array(total_ray_list,dtype=tuple).T.tolist()
        df = pd.DataFrame(total_ray_list,index=self.titles).T
        
        csv_output = 'Ray_Data.csv'
        with open(csv_output, 'w') as f:
            df.to_csv(f,index=False)