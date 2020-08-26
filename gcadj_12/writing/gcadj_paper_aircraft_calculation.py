import pandas as pd
import numpy as np


def Barometric_formula(P):
    Pb = 101325.0
    Tb = 288.15
    hb = 0
    R = 8.3144598
    Lb = -0.0065
    g0 = 9.80665
    M = 0.0289644
    h = (Tb / np.exp(np.log(P / Pb) * R * Lb / (g0 * M)) - Tb) / Lb + hb
    return h

class Aircraft():

    def __init__(self, aircraft_site_name,
                 lat, lon, pri_co2,
                 pos_co2, aircraft_co2,
                 height):
        self.site_name = aircraft_site_name
        self.lat = lat
        self.lon = lon
        self.pri_co2 = pri_co2
        self.pos_co2 = pos_co2
        self.aircraft_co2 = aircraft_co2
        self.height = height
        
    def bin_for_500m(self):
        height_list = []
        prior_co2_list = []
        posterior_co2_list = []
        aircraft_co2_list = []
        for site_num in range(self.site_name):
            bnum = 0

            tmp_height_list = []
            tmp_prior_co2_list = []
            tmp_posterior_co2_list = []
            tmp_aircraft_co2_list = []
            while bnum != self.height[site_num].shape[0]:
                bheight = self.height[site_num][bnum]
                enum  = np.where(self.height[site_num]<=bheight+500)[0][-1]+1
                nheight = np.mean(self.height[site_num][bnum:enum])
                npri_co2 = np.mean(self.pri_co2[site_num][bnum:enum])
                npos_co2 = np.mean(self.pos_co2[site_num][bnum:enum])
                naircraft_co2 = np.mean(self.aircraft_co2[site_num][bnum:enum])

                tmp_height_list.append(nheight)
                tmp_prior_co2_list.append(npri_co2)
                tmp_posterior_co2_list.append(npos_co2)
                tmp_aircraft_co2_list.append(naircraft_co2)
            height_list.append(np.array(tmp_height_list))
            prior_co2_list.append(np.array(tmp_prior_co2_list))
            posterior_co2_list.append(np.array(tmp_posterior_co2_list))
            aircraft_co2_list.append(np.array(tmp_aircraft_co2_list))
        self.height = height_list
        self.pri_co2 = prior_co2_list
        self.pos_co2 = posterior_co2_list
        self.aircraft_co2 = aircraft_co2_list
aircraft_pd = pd.read_excel('/home/huoxiao/n1s1/huoxiao_share/paper_data/gc_aircraft.xlsx', index_col=0)
aircraft_pd_num = int(aircraft_pd.shape[1]/9)


prior_co2_list = []
posterior_co2_list = []
aircraft_co2_list = []
height_list = []
for col in range(aircraft_pd_num):
    tmp_pressure = np.array(aircraft_pd.iloc[:, 9*col+5][~np.isnan(aircraft_pd.iloc[:, 9*col+5])])
    tmp_obspack_co2 = np.array(aircraft_pd.iloc[:, 9 * col + 6][~np.isnan(aircraft_pd.iloc[:, 9 * col + 6])])
    tmp_prior_co2 = np.array(aircraft_pd.iloc[:, 9*col+7][~np.isnan(aircraft_pd.iloc[:, 9*col+7])])
    tmp_posterior_co2 = np.array(aircraft_pd.iloc[:, 9 * col + 8][~np.isnan(aircraft_pd.iloc[:, 9 * col + 8])])
    tmp_height = Barometric_formula(tmp_pressure)
    print(tmp_height)
    prior_co2_list.append(tmp_prior_co2[tmp_height.argsort()])
    posterior_co2_list.append(tmp_posterior_co2[tmp_height.argsort()])
    aircraft_co2_list.append(tmp_obspack_co2[tmp_height.argsort()])
    tmp_height.sort()
    height_list.append(tmp_height)
