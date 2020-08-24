import xarray as xr
from netCDF4 import Dataset
import numpy as np
import warnings

class OCO2():

    def __init__(self, lat2D, lon2D):
        self.lat2D = lat2D
        self.lon2D = lon2D


    def search_xco2(self, llat, mlat,
                    llon, mlon, xco2):

        lat_2D_ht_llat = self.lat2D[np.where(self.lat2D >= llat)]
        lon_2D_ht_llat = self.lon2D[np.where(self.lat2D >= llat)]
        xco2_2D_ht_llat = xco2[np.where(self.lat2D >= llat)]

        lon_2D_llat_2_mlat = lon_2D_ht_llat[np.where(lat_2D_ht_llat < mlat)]
        xco2_2D_llat_2_mlat = xco2_2D_ht_llat[np.where(lat_2D_ht_llat < mlat)]

        # 90W TO 45W
        lon_2D_ht_llon = lon_2D_llat_2_mlat[np.where(lon_2D_llat_2_mlat > llon)]
        xco2_2D_ht_llon = xco2_2D_llat_2_mlat[np.where(lon_2D_llat_2_mlat > llon)]

        xco2_2D_90W_2_45W = xco2_2D_ht_llon[np.where(lon_2D_ht_llon < mlon)]

        return xco2_2D_90W_2_45W

warnings.filterwarnings('error')

# oco2_xr = xr.open_dataset('C:\\Users\\37014\\Desktop\\paper_data\\GEOSChem.XCO2onOCO2grid.nc4',
#                           decode_times=False,  decode_cf=False, use_cftime=False)
oco2_nc = Dataset('C:\\Users\\37014\\Desktop\\paper_data\\GEOSChem.XCO2onOCO2grid.nc4', 'r')
data_name_list = ['prior_xco2_nadir',
                  'posterior_xco2_nadir',
                  'oco2_xco2_nadir']

# oco2_lat = oco2_xr.coords['lat']
# oco2_lon = oco2_xr.coords['lon']
oco2_lat = oco2_nc.variables['lat'][:]
oco2_lon = oco2_nc.variables['lon'][:]
oco2_lon2D, oco2_lat2D = np.meshgrid(oco2_lon, oco2_lat)

oco2 = OCO2(oco2_lat2D, oco2_lon2D)
xco2_out_list = []
hxco2_out_list = []
for ith in range(3):
    # xco2 = oco2_xr.variables[data_name_list[ith]]
    xco2 = oco2_nc.variables[data_name_list[ith]][...]
    # 45N to 60 N
    rxco2 = oco2.search_xco2(45, 60, -180, 180, np.array(xco2[0, ...]))
    xco2_out_list.append(rxco2[np.where(rxco2!=0)])
    # 15S to 60S, 90W to 45W
    rxco2 = oco2.search_xco2(-60, -15, -90, -45, np.array(xco2[0, ...]))
    xco2_out_list.append(rxco2[np.where(rxco2!=0)])
    # 15S to 60S, 90W to 45W
    rxco2 = oco2.search_xco2(-60, -5, 100, 180, np.array(xco2[0, ...]))

    xco2_out_list.append(rxco2[np.where(rxco2!=0)])
    #higher
    # 20W to 10E
#mean
prior_xco2_mean = np.mean(np.array(xco2_out_list[0:9:3]))
try:
    print(prior_xco2_mean)
except:
    print('prior mean error')
    exit()
posterior_xco2_mean = np.mean(np.array(xco2_out_list[1:10:3]))
try:
    print(posterior_xco2_mean)
except:
    print('posterior mean error')
    exit()
oco2_xco2_mean = np.mean(np.array(xco2_out_list[2:11:3]))
try:
    print(oco2_xco2_mean)
except:
    print('oco2 mean error')
    exit()