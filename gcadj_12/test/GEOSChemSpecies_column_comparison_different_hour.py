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
import calendar

root_dir = '/home/huoxiao/n1s1/huoxiao_share/paper_data/GEOSChem.XCO2_diff_hour.nc4'
gc_xco2_nc = Dataset(root_dir, 'r')

year_mean0 = np.mean(gc_xco2_nc.variables['prior_xco20'][:], axis=0)
year_mean4 = np.mean(gc_xco2_nc.variables['prior_xco24'][:], axis=0)
year_mean3 = np.mean(gc_xco2_nc.variables['prior_xco23'][:], axis=0)
year_mean6 = np.mean(gc_xco2_nc.variables['prior_xco26'][:], axis=0)
year_mean12 = np.mean(gc_xco2_nc.variables['prior_xco212'][:], axis=0)

sum_days = 0

# MNL = MaxNLocator(200)
# lats = gc_xco2_nc.variables['lat'][:]
# lons = gc_xco2_nc.variables['lon'][:]
#
# out_pri_xco2 = year_mean0
# out_pos_xco2 = year_mean0
#
# x, y = np.meshgrid(lons, lats)
#
# fig = plt.figure(figsize=(16, 9))
# projection = ccrs.PlateCarree(central_longitude=0)
# axes_class = (GeoAxes, dict(map_projection=projection))
#
# prior_list = [year_mean0, year_mean3, year_mean6, year_mean12]
# for i in range(4):
#     out_pri_xco2 = prior_list[i]
#     axgr = AxesGrid(fig, 221+i, axes_class=axes_class,
#                     nrows_ncols=(1, 1),
#                     axes_pad=0.6,
#                     cbar_location='right',
#                     cbar_mode='each',
#                     cbar_pad=0.2,
#                     cbar_size='3%',
#                     label_mode='')  # note the empty label_mode
#     for i, ax in enumerate(axgr):
#         ax.set_extent([0,358,-90,91],crs=ccrs.PlateCarree())  # 120E ~30W
#         ax.set_xticks(np.arange(45,355,45),crs=ccrs.PlateCarree())
#         ax.set_yticks(np.arange(-90,120,30),crs=ccrs.PlateCarree())
#         ax.grid(linestyle = ":",color='black')
#
#         ax.add_feature(cfeature.BORDERS)
#         # map.add_feature(cfeature.OCEAN)
#         ax.add_feature(cfeature.COASTLINE)
#         p = ax.contourf(x, y, out_pri_xco2*10**6,
#                         transform=projection,
#                         cmap=get_cmap("jet"),
#                         locator=MNL,
#                         vmin=400, vmax=412)
#
#     cbar = axgr.cbar_axes[0].colorbar(p, ticks=np.arange(400, 413, 2))
#     # cbar.set_clim(400, 413)
#
#
#
# plt.show()
# sum_days = 0
# for month in range(12):
#     cur_days = calendar.monthrange(2018, month+1)[1]
#     year_mean0 = np.mean(gc_xco2_nc.variables['prior_xco20'][sum_days:sum_days+cur_days, ...], axis=0)
#     year_mean4 = np.mean(gc_xco2_nc.variables['prior_xco24'][sum_days:sum_days+cur_days, ...], axis=0)
#     year_mean3 = np.mean(gc_xco2_nc.variables['prior_xco23'][sum_days:sum_days+cur_days, ...], axis=0)
#     year_mean6 = np.mean(gc_xco2_nc.variables['prior_xco26'][sum_days:sum_days+cur_days, ...], axis=0)
#     year_mean12 = np.mean(gc_xco2_nc.variables['prior_xco212'][sum_days:sum_days+cur_days, ...], axis=0)
#
#     print('\n0-3 diff: ', (year_mean0*10**6-year_mean3*10**6).max(), (year_mean0*10**6-year_mean3*10**6).min(),
#           '0-6 diff: ', (year_mean0*10**6-year_mean6*10**6).max(), (year_mean0*10**6-year_mean6*10**6).min(),
#           '0-12 diff: ', (year_mean0*10**6-year_mean12*10**6).max(), (year_mean0*10**6-year_mean12*10**6).min())
#     sum_days += cur_days

#
# for week in range(50):
#
#     year_mean0 = np.mean(gc_xco2_nc.variables['prior_xco20'][week*7:week*7+7 ,...], axis=0)
#     year_mean4 = np.mean(gc_xco2_nc.variables['prior_xco24'][week*7:week*7+7, ...], axis=0)
#     year_mean3 = np.mean(gc_xco2_nc.variables['prior_xco23'][week*7:week*7+7, ...], axis=0)
#     year_mean6 = np.mean(gc_xco2_nc.variables['prior_xco26'][week*7:week*7+7, ...], axis=0)
#     year_mean12 = np.mean(gc_xco2_nc.variables['prior_xco212'][week*7:week*7+7, ...], axis=0)
#
#     print('\n0-3 diff: ', (year_mean0*10**6-year_mean3*10**6).max(), (year_mean0*10**6-year_mean3*10**6).min(),
#           '0-6 diff: ', (year_mean0*10**6-year_mean6*10**6).max(), (year_mean0*10**6-year_mean6*10**6).min(),
#           '0-12 diff: ', (year_mean0*10**6-year_mean12*10**6).max(), (year_mean0*10**6-year_mean12*10**6).min())
for day in range(365):

    year_mean0 = np.mean(gc_xco2_nc.variables['prior_xco20'][day:day+1 ,...], axis=0)
    year_mean4 = np.mean(gc_xco2_nc.variables['prior_xco24'][day:day+1, ...], axis=0)
    year_mean3 = np.mean(gc_xco2_nc.variables['prior_xco23'][day:day+1, ...], axis=0)
    year_mean6 = np.mean(gc_xco2_nc.variables['prior_xco26'][day:day+1, ...], axis=0)
    year_mean12 = np.mean(gc_xco2_nc.variables['prior_xco212'][day:day+1, ...], axis=0)

    print('\n0-3 diff: ', (year_mean0*10**6-year_mean3*10**6).max(), (year_mean0*10**6-year_mean3*10**6).min(),
          '0-6 diff: ', (year_mean0*10**6-year_mean6*10**6).max(), (year_mean0*10**6-year_mean6*10**6).min(),
          '0-12 diff: ', (year_mean0*10**6-year_mean12*10**6).max(), (year_mean0*10**6-year_mean12*10**6).min())