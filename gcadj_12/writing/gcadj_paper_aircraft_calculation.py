import pandas as pd
from netCDF4 import Dataset
import datetime
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.cm import get_cmap
from matplotlib.ticker import MaxNLocator


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
        for site_num in range(len(self.site_name)):

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
                bnum = enum
            height_list.append(np.array(tmp_height_list))
            prior_co2_list.append(np.array(tmp_prior_co2_list))
            posterior_co2_list.append(np.array(tmp_posterior_co2_list))
            aircraft_co2_list.append(np.array(tmp_aircraft_co2_list))
        self.height = height_list
        self.pri_co2 = prior_co2_list
        self.pos_co2 = posterior_co2_list
        self.aircraft_co2 = aircraft_co2_list

    def aircraft_plot(self, aircraft_index):

        fig, axes = plt.subplots(3, 4, figsize=(32, 16))
        for row in range(3):
            for col in range(4):
                ax = axes[row, col]
                print(aircraft_index[row+col], len(self.pri_co2))
                ax.plot(self.pri_co2[aircraft_index[row*4+col]], self.height[aircraft_index[row*4+col]])
                ax.plot(self.pos_co2[aircraft_index[row * 4 + col]], self.height[aircraft_index[row * 4 + col]])
                ax.plot(self.aircraft_co2[aircraft_index[row * 4 + col]], self.height[aircraft_index[row * 4 + col]])

    def aircraft_location_plot(self, aircraft_site_index):
        print(self.height[0])
        print(self.pri_co2[0])
        fig = plt.figure(figsize=(16, 9))
        projection = ccrs.PlateCarree(central_longitude=0)
        axes_class = (GeoAxes, dict(map_projection=projection))

        axgr = AxesGrid(fig, 111, axes_class=axes_class,
                        nrows_ncols=(1, 1),
                        axes_pad=0.6,
                        cbar_location='right',
                        cbar_mode=None,
                        cbar_pad=0.2,
                        cbar_size='3%',
                        label_mode='')  # note the empty label_mode
        for i, ax in enumerate(axgr):
            ax.set_extent([0, 358, -90, 91], crs=ccrs.PlateCarree())  # 120E ~30W
            ax.set_xticks(np.arange(45, 355, 45), crs=ccrs.PlateCarree())
            ax.set_yticks(np.arange(-90, 120, 30), crs=ccrs.PlateCarree())
            ax.grid(linestyle=":", color='black')
            ax.add_feature(cfeature.BORDERS)
            # map.add_feature(cfeature.OCEAN)
            ax.add_feature(cfeature.COASTLINE)


            ax.set_global()

            pri_bias_list_ht_2000 = []
            pri_bias_list_lt_2000 = []
            pos_bias_list_ht_2000 = []
            pos_bias_list_lt_2000 = []
            site_name_list = []
            for i in range(len(aircraft_site_index)):
                ax.text(lon_list[aircraft_site_index[i]],
                        lat_list[aircraft_site_index[i]],
                        aircraft_site_name[aircraft_site_index[i]][0:4])
                ax.scatter(self.lon[aircraft_site_index[i]], self.lat[aircraft_site_index[i]])

                #200
                pri_lt2000 = self.pri_co2[aircraft_site_index[i]][np.where(self.height[aircraft_site_index[i]]<2000)]
                pri_ht2000 = self.pri_co2[aircraft_site_index[i]][np.where(self.height[aircraft_site_index[i]]>=2000)]
                pos_lt2000 = self.pos_co2[aircraft_site_index[i]][np.where(self.height[aircraft_site_index[i]]<2000)]
                pos_ht2000 = self.pos_co2[aircraft_site_index[i]][np.where(self.height[aircraft_site_index[i]]>=2000)]
                aircraft_lt2000 = self.aircraft_co2[aircraft_site_index[i]][np.where(self.height[aircraft_site_index[i]]<2000)]
                aircraft_ht2000 = self.aircraft_co2[aircraft_site_index[i]][np.where(self.height[aircraft_site_index[i]]>=2000)]

                pri_bias_list_ht_2000.append(np.mean(pri_ht2000-aircraft_ht2000))
                pri_bias_list_lt_2000.append(np.mean(pri_lt2000 - aircraft_lt2000))
                pos_bias_list_ht_2000.append(np.mean(pos_ht2000-aircraft_ht2000))
                pos_bias_list_lt_2000.append(np.mean(pos_lt2000 - aircraft_lt2000))
                site_name_list.append(aircraft_site_name[aircraft_site_index[i]])
                # print(aircraft_site_name[aircraft_site_index[i]] + \
                # ' bias ht 2000', np.mean(pri_ht2000-aircraft_ht2000), np.mean(pos_ht2000-aircraft_ht2000))
                # print(aircraft_site_name[aircraft_site_index[i]] + \
                # ' bias lt 2000', np.mean(pri_lt2000 - aircraft_lt2000), np.mean(pos_lt2000 - aircraft_lt2000))
            site_name_array = np.array(site_name_list)[np.array(pri_bias_list_ht_2000).argsort()]
            pri_ht2000_array = np.array(pri_bias_list_ht_2000)
            pri_ht2000_array.sort()

            print(site_name_array)
            print(pri_ht2000_array)

            site_name_array = np.array(site_name_list)[np.array(pos_bias_list_ht_2000).argsort()]
            pos_ht2000_array = np.array(pos_bias_list_ht_2000)
            pos_ht2000_array.sort()
            print(site_name_array)
            print(pos_ht2000_array)
            # print('bias ht 2000', np.mean(np.array(pri_bias_list_ht_2000)), np.mean(np.array(pri_bias_list_lt_2000)),
            #       np.mean(np.array(pos_bias_list_ht_2000)), np.mean(np.array(pos_bias_list_lt_2000)))
            #
            # print('bias max min ', np.array(pri_bias_list_ht_2000).max(), np.array(pri_bias_list_ht_2000).min(),
            #       np.array(pos_bias_list_ht_2000).max(), np.array(pos_bias_list_ht_2000).min())
aircraft_pd = pd.read_excel('/home/huoxiao/n1s1/huoxiao_share/paper_data/gc_aircraft_obs_flag_gc.xlsx', index_col=0)
aircraft_pd_num = int(aircraft_pd.shape[1]/10)


prior_co2_list = []
posterior_co2_list = []
aircraft_co2_list = []
height_list = []
lat_list = []
lon_list = []
time_list = []
blh_list = []

aircraft_site_index = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
aircraft_site_name = ['Briggsdale', 'Offshore Cape May', 'Estevan Point', 'East Trout Lake',
                      'Homer', 'Park Falls', 'Offshore Portsmouth', 'Poker Flat',
                      'Rarotonga', 'Offshore Charleston', 'Southern Great Plains',
                      'Trinidad Head, California', 'West Branch']
for col in range(aircraft_pd_num):
    tmp_pressure = np.array(aircraft_pd.iloc[:, 10*col+5][~np.isnan(aircraft_pd.iloc[:, 10*col+5])])
    tmp_obspack_co2 = np.array(aircraft_pd.iloc[:, 10 * col + 6][~np.isnan(aircraft_pd.iloc[:, 10 * col + 6])])*10**6
    tmp_prior_co2 = np.array(aircraft_pd.iloc[:, 10*col+7][~np.isnan(aircraft_pd.iloc[:, 10*col+7])])*10**6
    tmp_posterior_co2 = np.array(aircraft_pd.iloc[:, 10 * col + 8][~np.isnan(aircraft_pd.iloc[:, 10 * col + 8])])*10**6
    lat_list.append(np.mean(aircraft_pd.iloc[:, 10 * col + 2][~np.isnan(aircraft_pd.iloc[:, 10 * col + 2])]))
    lon_list.append(np.mean(aircraft_pd.iloc[:, 10 * col + 1][~np.isnan(aircraft_pd.iloc[:, 10 * col + 1])]))
    time_list.append(aircraft_pd.iloc[:, 10*col+0][~np.isnan(aircraft_pd.iloc[:, 10*col+5])])
    blh_list.append(np.mean(aircraft_pd.iloc[:, 10 * col + 9][~np.isnan(aircraft_pd.iloc[:, 10 * col + 1])]))

    tmp_height = Barometric_formula(tmp_pressure)
    prior_co2_list.append(tmp_prior_co2[tmp_height.argsort()])
    posterior_co2_list.append(tmp_posterior_co2[tmp_height.argsort()])
    aircraft_co2_list.append(tmp_obspack_co2[tmp_height.argsort()])
    tmp_height.sort()
    height_list.append(tmp_height)

CAircraft = Aircraft(aircraft_site_name, lat_list,
                     lon_list, prior_co2_list,
                     posterior_co2_list, aircraft_co2_list,
                     height_list)
CAircraft.bin_for_500m()
# CAircraft.aircraft_plot(aircraft_site_index)
CAircraft.aircraft_location_plot(aircraft_site_index)
print(blh_list)

plt.show()

#[2013.3932543122366, 472.0210489069374, 645.2858927572927, 947.2610752146732, 1119.05687445925, 997.1165771484375, 789.0087252862439, 697.1826794757935, 610.6280536651611, 671.3800738609001, 1170.2290109773962, 442.62534791231155, 1135.7017129968715]