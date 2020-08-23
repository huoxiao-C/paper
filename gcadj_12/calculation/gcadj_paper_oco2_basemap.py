
from netCDF4 import Dataset, num2date, date2num
import datetime
import subprocess
import numpy as np

import time

import ESMF


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


def get_ij(lon, lat, delta_x, delta_y, NX, NY):

    ''''
    '''
    tlon = int((lon+180.0)/delta_x+1.5)

    tlat = int((lat + 90.0) / delta_y + 1.5)

    if( tlon > NX ):
        tlon = tlon-NX

    return tlon, tlat

def pressure_weigt_function(pres, surfp, level, lat_num, lon_num):
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

def gc_interp_tccon(gc_co2_native, gc_level, gc_pres, tccon_level, tccon_pres, tccon_surfp):
    ''''
    # map geos-chem profile concentration to tccon profile concentration
    '''
    hinterpz = np.zeros((gc_level, tccon_level))
    for Ltccon in range(tccon_level):
        for Lgc in range(gc_level-1):
            low = gc_pres[Lgc+1]
            hi = gc_pres[Lgc]
            if(tccon_pres[Ltccon]<=hi and tccon_pres[Ltccon]>low):
                diff = hi-low
                hinterpz[Lgc+1, Ltccon] = (hi-tccon_pres[Ltccon])/diff
                hinterpz[Lgc, Ltccon] = (tccon_pres[Ltccon]-low)/diff

    for Ltccon in range(tccon_level):
        if(tccon_pres[Ltccon]>gc_pres[0]):
            hinterpz[0, Ltccon] = 1
            hinterpz[1:gc_level, Ltccon] = 0

    gc_co2 = np.dot(gc_co2_native, hinterpz)
    return gc_co2


def gc_co2_2_xco2(gc_co2, gc_level, gc_pres, oco2_nc_data, oco2_level, dataith):
    '''

    :param gc_co2:
    :param gc_level:
    :param gc_pres:
    :param oco2_nc_data:
    :return:
    '''
    oco2_pres = oco2_nc_data.variables['pressure'][dataith, :]
    oco2_surp_pres = oco2_nc_data.variables['surface_pressure'][dataith]
    gx_co2_inter_oco2_level = gc_interp_tccon(gc_co2, gc_level, gc_pres, oco2_level, oco2_pres, oco2_surp_pres)
    gc_xco2 = oco2_nc_data.variables['xco2_apriori'][dataith]+np.sum((gx_co2_inter_oco2_level-oco2_nc_data.variables['co2_profile_apriori'][dataith, :])*oco2_nc_data.variables['xco2_avg_kernel'][dataith, :]\
              *oco2_nc_data.variables['xco2_pressure_weighting_function'][dataith, :])
    return gc_xco2

begin_time_str = "20180101000000"
end_time_str = "20190101000000"
begin_time = datetime.datetime.strptime(begin_time_str, '%Y%m%d%H%M%S')
end_time = datetime.datetime.strptime(end_time_str, '%Y%m%d%H%M%S')

gc_pri_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/'
gc_pos_dir = '/share/nas1_share1/huoxiao_share/'
ct_dir = '/share/nas1_share1/huoxiao_share/CO2_data/CarbonTracker/CO2_Concentration/'
oco2_dir = '/share/nas1_share1/huoxiao_share/CO2_data/Observation_data/satellite_on_gc_2x25grid/'

oco2_filename_list = ['OCO2_CO2_Litev9_ND_Observation_Error_Covariance', 'OCO2_CO2_Litev9_ND_plus_GlintLand', 'OCO2_CO2_Litev9_ND_plus_GlintLand_plus_GlintOcean']
gc_filename_list = ['GEOSChem_Assimilation_Litev9_Nadir_Obsm_v2', 'GEOSChem_Assimilation_Litev9_Nadir_plus_GlintLand_Obs', \
                    'GEOSChem_Assimilation_Litev9_Nadir_plus_GlintLand_plus_GlintOcean_Obs']
oco2_filename_bre = ['nadir', 'nadir_plus_glintland', 'nd_plus_glintland_plus_glintocean']
#sinal whether have different assimilation data
NADIR = 1
NADIR_PLUS_GLINTLAND = 0
NADIR_PLUS_GLINTLAND_GLINTOCEAN = 0
oco2_obs_flag_array = np.array([NADIR, NADIR_PLUS_GLINTLAND, NADIR_PLUS_GLINTLAND_GLINTOCEAN])
#initialization
delta_x = 2.5
delta_y = 2
NX = 144
NY = 91
gc_level = 47
oco2_level = 20
lat_num = 91
lon_num = 144
ct_level = 25
ct_lat_num = 90
ct_lon_num =120
#array
prior_co2_array = np.zeros((3, lat_num, lon_num))
posterior_gc_co2_array = np.zeros((3, lat_num, lon_num))
oco2_co2_array = np.zeros((3, lat_num, lon_num))
ct_co2_array = np.zeros((3, lat_num, lon_num))
co2_num_array = np.zeros((3, lat_num, lon_num))

#
pbtime = time.time()

#mean difference between model observation and satellite observation

while( end_time != begin_time ):
    print('stage1', datetime.datetime.strftime(begin_time, '%Y%m%d%H%M%S'))

    for oco2_obsith in range(oco2_obs_flag_array.shape[0]):
        if oco2_obs_flag_array[oco2_obsith] == 0:
            continue
        oco2_file_path = oco2_dir+oco2_filename_list[oco2_obsith]+'/oco2_'+ begin_time.strftime('%Y%m%d%H%M%S') + '.nc'

        if not subprocess.call('test -f '+oco2_file_path, shell=True):

            oco2_nc_data = Dataset(oco2_file_path, 'r')

            tmp_xco2_file = oco2_nc_data.variables['xco2'][:]
            tmp_lat_file  = oco2_nc_data.variables['latitude'][:]
            tmp_lon_file  = oco2_nc_data.variables['longitude'][:]
            tmp_xco2_uncert_file = oco2_nc_data.variables['xco2_uncert'][:]


            #prior

            gc_met_filename = 'GEOSChem.StateMet.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                            .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour, begin_time.minute)

            gc_pri_spc_filename = 'GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                            .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour, begin_time.minute)


            ct_filename = '{:0>4d}/CT2019.molefrac_glb3x2_{:0>4d}-{:0>2d}-{:0>2d}.nc'. \
                format(begin_time.year, begin_time.year, begin_time.month, begin_time.day)

            gc_pri_spc_path = gc_pri_dir+gc_pri_spc_filename
            gc_met_path = gc_pri_dir+gc_met_filename
            ct_file_path = ct_dir+ct_filename

            gc_pri_spc_nc_data = Dataset(gc_pri_spc_path, 'r')
            gc_met_nc_data = Dataset(gc_met_path, 'r')
            ct_nc_data = Dataset(ct_file_path, 'r')

            ct_lons = ct_nc_data.variables['longitude'][:]
            ct_lats = ct_nc_data.variables['latitude'][:]
            gc_lons = gc_met_nc_data.variables['lon'][:]
            gc_lats = gc_met_nc_data.variables['lat'][:]
            ct_pre = ct_nc_data.variables['pressure'][int(begin_time.hour/3), :, :, :]
            ct_surfp = ct_nc_data.variables['pressure'][int(begin_time.hour/3), 0, :, :]
            ct_h = pressure_weigt_function(ct_pre, ct_surfp, ct_level, ct_lat_num, ct_lon_num)
            ct_xco2 = np.sum(ct_h[:, :, :]*ct_nc_data.variables['co2'][int(begin_time.hour/3), :, :, :], axis=0)
            ct2gc_xco2 = ct_concentration2gc_concentration(ct_xco2, ct_lons, ct_lats, gc_lons, gc_lats)

            h  = pressure_weigt_function(gc_met_nc_data.variables['Met_PMID'][0, :, :, :], \
                                        gc_met_nc_data.variables['Met_PMID'][0, 0, :, :], gc_level, NY, NX)

            gc_pri_co2 = gc_pri_spc_nc_data.variables['SpeciesConc_CO2'][0, :, :, :]
            #posterior

            gc_pos_spc_filename = '/GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                            .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour, begin_time.minute)

            gc_pos_spc_path = gc_pos_dir+gc_filename_list[oco2_obsith]+gc_pos_spc_filename



            gc_pos_spc_nc_data = Dataset(gc_pos_spc_path, 'r')
            gc_pos_co2 = gc_pos_spc_nc_data.variables['SpeciesConc_CO2'][0, :, :, :]

            for oco2ith in range(tmp_xco2_file.shape[0]):
                lon_index, lat_index = get_ij(tmp_lon_file[oco2ith], tmp_lat_file[oco2ith], delta_x, delta_y, NX, NY)
                lon_index = lon_index-1
                lat_index = lat_index-1
                gc_pri_xco2 =  gc_co2_2_xco2(gc_pri_co2[:, lat_index, lon_index], gc_level, \
                                             gc_met_nc_data.variables['Met_PMID'][0, :, lat_index, lon_index]*100, \
                                             oco2_nc_data, oco2_level, oco2ith)
                gc_pos_xco2 =  gc_co2_2_xco2(gc_pos_co2[:, lat_index, lon_index], gc_level, \
                                             gc_met_nc_data.variables['Met_PMID'][0, :, lat_index, lon_index]*100, \
                                             oco2_nc_data, oco2_level, oco2ith)
                prior_co2_array[oco2_obsith, lat_index, lon_index] += gc_pri_xco2
                oco2_co2_array[oco2_obsith, lat_index, lon_index] += tmp_xco2_file[oco2ith]
                posterior_gc_co2_array[oco2_obsith, lat_index, lon_index] += gc_pos_xco2
                ct_co2_array[oco2_obsith, lat_index, lon_index] += ct2gc_xco2[lat_index, lon_index]
                co2_num_array[oco2_obsith, lat_index, lon_index] += 1


            #close file
            gc_pri_spc_nc_data.close()
            gc_met_nc_data.close()
            gc_pos_spc_nc_data.close()
            ct_nc_data.close()
            oco2_nc_data.close()

    begin_time = begin_time + datetime.timedelta(minutes=60)

#average
prior_co2_array[np.where(co2_num_array!=0)] /= co2_num_array[np.where(co2_num_array!=0)]
oco2_co2_array[np.where(co2_num_array!=0)] /= co2_num_array[np.where(co2_num_array!=0)]
posterior_gc_co2_array[np.where(co2_num_array!=0)] /= co2_num_array[np.where(co2_num_array!=0)]
ct_co2_array[np.where(co2_num_array!=0)] /= co2_num_array[np.where(co2_num_array!=0)]


out_file_path = '/home/huoxiao/paper_data/GEOSChem.XCO2onOCO2grid.nc4'
out_oco2 = Dataset(out_file_path, 'w')
Dataset.createDimension(out_oco2, dimname='time', size=None)
Dataset.createDimension(out_oco2, dimname='lev', size=47)
Dataset.createDimension(out_oco2, dimname='lat', size=lat_num)
Dataset.createDimension(out_oco2, dimname='lon', size=lon_num)
Dataset.createDimension(out_oco2, dimname='DateStrLen', size=23)

Times_out = Dataset.createVariable(out_oco2, 'Times', datatype=np.float64, dimensions=('time'))
Times_out.units = "seconds since 1970-01-01 00:00:00"
Times_out[:] = date2num(begin_time, "seconds since 1970-01-01 00:00:00")

level_out = Dataset.createVariable(out_oco2, 'level', datatype=np.float32, dimensions=('lev'))
level_out[:] = np.arange(1, 48, 1)

lat_out = Dataset.createVariable(out_oco2, 'lat', datatype=np.float32, dimensions=('lat'))
lat_out[:] = gc_lats

lon_out = Dataset.createVariable(out_oco2, 'lon', datatype=np.float32, dimensions=('lon'))
lon_out[:] = gc_lons

for oco2_obsith in range(oco2_obs_flag_array.shape[0]):
    if oco2_obs_flag_array[oco2_obsith] == 0:
        continue
    oco2_xco2_out = Dataset.createVariable(out_oco2, 'oco2_xco2_'+oco2_filename_bre[oco2_obsith], \
                                           datatype=np.float32, dimensions=('time', 'lat', 'lon'))
    oco2_xco2_out[...] = oco2_co2_array[oco2_obsith, ...]


    co2_apri_out = Dataset.createVariable(out_oco2, 'prior_xco2_'+oco2_filename_bre[oco2_obsith], \
                                          datatype=np.float32,dimensions=('time','lat', 'lon'))
    co2_apri_out[...] = prior_co2_array[oco2_obsith, ...]

    co2_pos_out = Dataset.createVariable(out_oco2, 'posterior_xco2_'+oco2_filename_bre[oco2_obsith], \
                                          datatype=np.float32,dimensions=('time', 'lat', 'lon'))
    co2_pos_out[...] = posterior_gc_co2_array[oco2_obsith, ...]

    ct_pos_out = Dataset.createVariable(out_oco2, 'ct_xco2_'+oco2_filename_bre[oco2_obsith], \
                                          datatype=np.float32,dimensions=('time', 'lat', 'lon'))
    ct_pos_out[...] = ct_co2_array[oco2_obsith, ...]
