
from netCDF4 import Dataset, num2date
import datetime
import numpy as np
import ESMF

def files_string2list(files_string):
    '''
    将字符串变为分割为文件名
    :param files_string:
    :return:文件名列表返回
    '''
    files_list = []
    temp = ''
    for i in range(len(files_string)):
        if files_string[i] != '\n':
            temp += files_string[i]
        else:
            files_list.append(temp)
            temp = ''
    return files_list

def pressure_weigt_function(gc_pre, q, Mdry, g, level, surface_pre):
    ''''

    '''

    # cbar = (1-q[1:level])/(g*Mdry)
    # gc_pre[0] = surface_pre
    # deltap = gc_pre[0:level-1]-gc_pre[1:level]
    # deltap[0] = gc_pre[1]-surface_pre
    # hdot = cbar*deltap/(np.sum(deltap*cbar))
    #
    # fi = 1.0/2.0
    # fs = (surface_pre-gc_pre[0])/(gc_pre[0]-gc_pre[1])
    #
    # h = np.zeros(level)
    #
    # for i in range(level):
    #     if i == 0:
    #         h[i] = fs*fi*hdot[i]
    #     elif i == 1:
    #         h[i] = fi*hdot[i]+(1-fs)*hdot[i-1]
    #     elif i == level-1:
    #         h[i] = (1-fi)*hdot[i-1]
    #     else:
    #         h[i] = fi*hdot[i]+(1-fi)*hdot[i-1]
    # return h

    cbar = (1-q)/(g*Mdry)
    gc_pre[0, :, :] = surface_pre[:, :]
    deltap = gc_pre[0:level, :, :]-gc_pre[1:level+1, :, :]
    h = cbar*deltap/(np.sum(deltap*cbar, axis=0))
    return h


def ct_pressure_weigt_function(pres, surfp, level, lat_num, lon_num):
    ''''

    '''
    h = np.zeros((level, lat_num, lon_num))
    for i in range(level):
        if i == 0:
            h[i, :, :] = np.abs(-pres[i, :, :] + (pres[i + 1, :, :] - pres[i, :, :]) / np.log(pres[i + 1, :, :] / pres[i, :, :]))/surfp[:, :]
        elif i == level-1:
            h[i, :, :] = np.abs(pres[i, :, :]-(pres[i, :, :] - pres[i - 1, :, :]) / np.log(pres[i, :, :] / pres[i - 1, :, :]))/surfp[:, :]
        else:
            h[i, :, :] = np.abs(-pres[i, :, :] + (pres[i + 1, :, :] - pres[i, :, :]) / np.log(pres[i + 1, :, :] / \
                         pres[i, :, :])+pres[i, :, :]-(pres[i, :, :] - pres[i - 1, :, :]) / np.log(pres[i, :, :] / pres[i - 1, :, :]))/surfp[:, :]

    return h

def ct_concentration2gc_concentration(ct_xco2, ct_lon, ct_lat, gc_lon, gc_lat):
    '''

    :param ct_xco2:
    :param ct_lon:
    :param ct_lat:
    :param gc_lon:
    :param gc_lat:
    :return:
    '''
    ct_lon_shape = ct_lon.shape[0]
    ct_lat_shape = ct_lat.shape[0]
    gc_lon_shape = gc_lon.shape[0]
    gc_lat_shape = gc_lat.shape[0]

    sourcegrid = ESMF.Grid(np.array([ct_lon_shape, ct_lat_shape]), staggerloc=ESMF.StaggerLoc.CENTER,
                           coord_sys=ESMF.CoordSys.SPH_DEG)
    destgrid = ESMF.Grid(np.array([gc_lon_shape, gc_lat_shape]), staggerloc=ESMF.StaggerLoc.CENTER,
                         coord_sys=ESMF.CoordSys.SPH_DEG)

    source_lon = sourcegrid.get_coords(0)
    source_lat = sourcegrid.get_coords(1)

    dest_lon = destgrid.get_coords(0)
    dest_lat = destgrid.get_coords(1)

    # transpose　for broadcast
    source_lon = source_lon.transpose()
    dest_lon = dest_lon.transpose()
    source_lon[...] = ct_lon
    source_lat[...] = ct_lat
    dest_lon[...] = gc_lon
    dest_lat[...] = gc_lat
    # reset
    source_lon = source_lon.transpose()
    dest_lon = dest_lon.transpose()

    sourcefield = ESMF.Field(sourcegrid, name='carbontracker co2')
    destfield = ESMF.Field(destgrid, name='geos-chem co2')

    sourcefield.data[...] = ct_xco2.transpose()
    regrid_blr = ESMF.Regrid(sourcefield, destfield, regrid_method=ESMF.RegridMethod.BILINEAR,
                             unmapped_action=ESMF.UnmappedAction.IGNORE)
    regrid_nst = ESMF.Regrid(sourcefield, destfield, regrid_method=ESMF.RegridMethod.NEAREST_STOD,
                             unmapped_action=ESMF.UnmappedAction.IGNORE)
    destfield_blr = regrid_blr(sourcefield, destfield)
    destfield_nst = regrid_nst(sourcefield, destfield)

    destfield_blr.data[np.where(destfield_blr.data == 0)] = destfield_nst.data[np.where(destfield_blr.data == 0)]
    ct_hor_xco2 = destfield_blr.data.transpose()
    return ct_hor_xco2

#initialization
root_dir = '/share/nas1_share1'
begin_time_str = '20180101000000'
end_time_str = "20190101000000"
begin_time = datetime.datetime.strptime(begin_time_str, '%Y%m%d%H%M%S')
ct_begin_time_date = datetime.datetime.strptime(begin_time_str, '%Y%m%d%H%M%S')
end_time = datetime.datetime.strptime(end_time_str, '%Y%m%d%H%M%S')
ct_end_time_date = datetime.datetime.strptime(end_time_str, '%Y%m%d%H%M%S')

days = (end_time-begin_time).days
lons      = np.zeros((91, 144))
lats      = np.zeros((91, 144))
xco2_ori0      = np.zeros((91, 144), dtype=float)
xco2_ass0      = np.zeros((91, 144), dtype=float)
xco2_diff0      = np.zeros((91, 144), dtype=float)

xco2_ori3      = np.zeros((91, 144), dtype=float)
xco2_ass3      = np.zeros((91, 144), dtype=float)
xco2_diff3      = np.zeros((91, 144), dtype=float)

xco2_ori4      = np.zeros((91, 144), dtype=float)
xco2_ass4      = np.zeros((91, 144), dtype=float)
xco2_diff4      = np.zeros((91, 144), dtype=float)

xco2_ori6      = np.zeros((91, 144), dtype=float)
xco2_ass6      = np.zeros((91, 144), dtype=float)
xco2_diff6      = np.zeros((91, 144), dtype=float)

xco2_ori12      = np.zeros((91, 144), dtype=float)
xco2_ass12      = np.zeros((91, 144), dtype=float)
xco2_diff12      = np.zeros((91, 144), dtype=float)
#gc
Mdry = 29.0
level = 47
g = 9.8
tot_num = 0

xco2_ori_day_mean_list = []
xco2_ass_day_mean_list = []
xco2_diff_day_mean_list = []
#calendar.monthrange(2016, month)[1]+1

year_xco2_pri0 = np.zeros((days, 91, 144))
year_xco2_ass0 = np.zeros((days, 91, 144))

year_xco2_pri3 = np.zeros((days, 91, 144))
year_xco2_ass3 = np.zeros((days, 91, 144))

year_xco2_pri4 = np.zeros((days, 91, 144))
year_xco2_ass4 = np.zeros((days, 91, 144))

year_xco2_pri6 = np.zeros((days, 91, 144))
year_xco2_ass6 = np.zeros((days, 91, 144))

year_xco2_pri12 = np.zeros((days, 91, 144))
year_xco2_ass12 = np.zeros((days, 91, 144))
for day in range(days):
    for hour in range(24):
        print(datetime.datetime.strftime(begin_time, '%Y-%m-%d %H:%M:%S'))

        gcrspc_bac_dir = root_dir+'/huoxiao_share/GEOSChem_Simulation_test/GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4'\
                        .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour)
        gcrspc_dir = root_dir+'/huoxiao_share/GEOSChem_Assimilation_Litev9_Nadir_Obsm_v2/GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4'\
            .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour)

        gcmet_dir = root_dir+'/huoxiao_share/GEOSChem_Simulation_test/GEOSChem.StateMet.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4'\
            .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour)
        print(gcmet_dir)

        gcmet = Dataset(gcmet_dir, 'r')
        gcspc = Dataset(gcrspc_dir, 'r')
        gcspc_bac = Dataset(gcrspc_bac_dir, 'r')

        lat_num =91
        lon_num = 144
        ver_hy = gcmet.variables['hybi'][:].shape[0]
        gc_pre = (gcmet.variables['hyai'][:] + np.kron(gcmet.variables['Met_PSC2WET'][0, :, :], gcmet.variables['hybi'][:]).reshape(lat_num, lon_num, ver_hy)).transpose(2, 0, 1)
        surface_pre = gcmet.variables['Met_PSC2WET'][0, :, :]

        q = gcmet.variables['Met_SPHU'][0, :, :, :] * 1e-3
        h = pressure_weigt_function(gc_pre, q, Mdry, g, level, surface_pre)

        xco2_ass0 += np.sum(h[:, :, :]*(gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)

        xco2_ori0 += np.sum(h[:, :, :]*(gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
        xco2_diff0 += np.sum(h[:, :, :]*(gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)-np.sum(h[:, :, :]*(gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)

        if hour%4==0:
            xco2_ass4 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)

            xco2_ori4 += np.sum(h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
            xco2_diff4 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0) - np.sum(
                h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
        if hour%3==0:
            xco2_ass3 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)

            xco2_ori3 += np.sum(h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
            xco2_diff3 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0) - np.sum(
                h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
        if hour%6==0:
            xco2_ass6 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)

            xco2_ori6 += np.sum(h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
            xco2_diff6 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0) - np.sum(
                h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
        if hour%12==0:
            xco2_ass12 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)

            xco2_ori12 += np.sum(h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
            xco2_diff12 += np.sum(h[:, :, :] * (gcspc.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0) - np.sum(
                h[:, :, :] * (gcspc_bac.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)

        tot_num += 1
        begin_time += datetime.timedelta(hours=1)

    year_xco2_pri0[day] = xco2_ori0/24
    year_xco2_ass0[day] = xco2_ass0/24
    year_xco2_pri4[day] = xco2_ori4/6
    year_xco2_ass4[day] = xco2_ass4/6
    year_xco2_pri3[day] = xco2_ori3/8
    year_xco2_ass3[day] = xco2_ass3/8
    year_xco2_pri6[day] = xco2_ori6/4
    year_xco2_ass6[day] = xco2_ass6/4
    year_xco2_pri12[day] = xco2_ori12/2
    year_xco2_ass12[day] = xco2_ass12/2
    #reset
    xco2_ass0 = 0
    xco2_diff0 = 0
    xco2_ori0 = 0

    xco2_ass4 = 0
    xco2_diff4 = 0
    xco2_ori4 = 0

    xco2_ass3 = 0
    xco2_diff3 = 0
    xco2_ori3 = 0

    xco2_ass6 = 0
    xco2_diff6 = 0
    xco2_ori6 = 0

    xco2_ass12 = 0
    xco2_diff12 = 0
    xco2_ori12 = 0
    #plot
    # xco2 = gc.variables['SpeciesRst_CO2'][0, 0, :, :]




gc_lons = gcmet.variables['lon'][:]
gc_lats = gcmet.variables['lat'][:]
#
# #convert catbontracker grid to geos-chem grid
ct_dir = root_dir+'/huoxiao_share/CO2_data/CarbonTracker/CO2_Concentration/'
ct_level = 25
ct_lat_num = 90
ct_lon_num =120

#method 2 vertical concentration vary linearly

ct_xco2 = np.zeros((ct_lat_num, ct_lon_num))
year_ct_xco2 = np.zeros((days, 91, 144))
for ct_day in range(days):
    filename = '{:0>4d}/CT2019.molefrac_glb3x2_{:0>4d}-{:0>2d}-{:0>2d}.nc'. \
        format(ct_begin_time_date.year, ct_begin_time_date.year, ct_begin_time_date.month, ct_begin_time_date.day)

    ct_file_path = ct_dir + filename
    print(ct_file_path)
    ct_nc_data = Dataset(ct_file_path, 'r')
    for ct_num in range(8):
        ct_pre  = ct_nc_data.variables['pressure'][ct_num, :, :, :]
        ct_surfp = ct_nc_data.variables['pressure'][ct_num, 0, :, :]
        h = ct_pressure_weigt_function(ct_pre, ct_surfp, ct_level, ct_lat_num, ct_lon_num)

        ct_xco2 += np.sum(h[:, :, :]*(ct_nc_data.variables['co2'][ct_num, :, :, :]), axis=0)

        ct_begin_time_date += datetime.timedelta(hours=3)
    ct_xco2 = ct_xco2/8*10**-6
    ct_lons = ct_nc_data.variables['longitude'][:]
    ct_lats = ct_nc_data.variables['latitude'][:]
    ct2gc_xco2 = ct_concentration2gc_concentration(ct_xco2, ct_lons, ct_lats, gc_lons, gc_lats)
    year_ct_xco2[ct_day] = ct2gc_xco2
    ct_xco2 = 0.0

#output
out_file_path = '/home/huoxiao/paper_data/GEOSChem.XCO2_diff_hour.nc4'
print('output file: ', out_file_path)
out_oco2 = Dataset(out_file_path, 'w')
Dataset.createDimension(out_oco2, dimname='time', size=None)
Dataset.createDimension(out_oco2, dimname='lev', size=47)
Dataset.createDimension(out_oco2, dimname='lat', size=lat_num)
Dataset.createDimension(out_oco2, dimname='lon', size=lon_num)
Dataset.createDimension(out_oco2, dimname='DateStrLen', size=23)

Times_out = Dataset.createVariable(out_oco2, 'Times', datatype=np.float64, dimensions=('time'))
Times_out.units = "days since 2018-01-01"
Times_out[...] = np.arange(day)

lat_out = Dataset.createVariable(out_oco2, 'lat', datatype=np.float32, dimensions=('lat'))
lat_out[:] = gc_lats

lon_out = Dataset.createVariable(out_oco2, 'lon', datatype=np.float32, dimensions=('lon'))
lon_out[:] = gc_lons

#oco2
pri_xco2_out = Dataset.createVariable(out_oco2, 'prior_xco20', \
                                       datatype=np.float32, dimensions=('time', 'lat', 'lon'))
pri_xco2_out[...] = year_xco2_pri0



pos_xco2_out = Dataset.createVariable(out_oco2, 'posterior_xco20', \
                                      datatype=np.float32,dimensions=('time', 'lat', 'lon'))
pos_xco2_out[...] = year_xco2_ass0

#4
pri_xco2_out = Dataset.createVariable(out_oco2, 'prior_xco24', \
                                       datatype=np.float32, dimensions=('time', 'lat', 'lon'))
pri_xco2_out[...] = year_xco2_pri4

pos_xco2_out = Dataset.createVariable(out_oco2, 'posterior_xco24', \
                                      datatype=np.float32,dimensions=('time', 'lat', 'lon'))
pos_xco2_out[...] = year_xco2_ass4

#3
pri_xco2_out = Dataset.createVariable(out_oco2, 'prior_xco23', \
                                       datatype=np.float32, dimensions=('time', 'lat', 'lon'))
pri_xco2_out[...] = year_xco2_pri3



pos_xco2_out = Dataset.createVariable(out_oco2, 'posterior_xco23', \
                                      datatype=np.float32,dimensions=('time', 'lat', 'lon'))
pos_xco2_out[...] = year_xco2_ass3
#6
pri_xco2_out = Dataset.createVariable(out_oco2, 'prior_xco26', \
                                       datatype=np.float32, dimensions=('time', 'lat', 'lon'))
pri_xco2_out[...] = year_xco2_pri6



pos_xco2_out = Dataset.createVariable(out_oco2, 'posterior_xco26', \
                                      datatype=np.float32,dimensions=('time', 'lat', 'lon'))
pos_xco2_out[...] = year_xco2_ass6

#12
pri_xco2_out = Dataset.createVariable(out_oco2, 'prior_xco212', \
                                       datatype=np.float32, dimensions=('time', 'lat', 'lon'))
pri_xco2_out[...] = year_xco2_pri12



pos_xco2_out = Dataset.createVariable(out_oco2, 'posterior_xco212', \
                                      datatype=np.float32,dimensions=('time', 'lat', 'lon'))
pos_xco2_out[...] = year_xco2_ass12

ct_xco2_out = Dataset.createVariable(out_oco2, 'ct_xco2', \
                                      datatype=np.float32,dimensions=('time', 'lat', 'lon'))
ct_xco2_out[...] = year_ct_xco2
