import xarray as xr
from netCDF4 import Dataset
import numpy as np

root_dir = '/home/huoxiao/n1s1/huoxiao_share/GEOSChem.XCO2_diff_hour.nc4'
gc_xco2_nc = Dataset(root_dir, 'r')

year_mean0 = np.mean(gc_xco2_nc.variables['prior_xco20'][:], axis=0)
year_mean4 = np.mean(gc_xco2_nc.variables['prior_xco24'][:], axis=0)
year_mean3 = np.mean(gc_xco2_nc.variables['prior_xco23'][:], axis=0)
year_mean6 = np.mean(gc_xco2_nc.variables['prior_xco26'][:], axis=0)
year_mean12 = np.mean(gc_xco2_nc.variables['prior_xco26'][:], axis=0)