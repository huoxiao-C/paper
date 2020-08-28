from netCDF4 import Dataset, date2num, num2date
import os, time
import numpy as np
import datetime
from fnmatch import fnmatch
import pandas as pd

def files_string2list(files_string):
    '''
    灏嗗瓧绗︿覆鍙樹负strptime鍒嗗壊涓烘枃浠跺悕
    :param files_string:
    :return:鏂囦欢鍚嶅垪琛ㄨ繑鍥?    '''
    files_list = []
    temp = ''
    for i in range(len(files_string)):
        if files_string[i] != '\n':

            temp += files_string[i]
        else:
            files_list.append(temp)
            temp = ''
    return files_list

def get_ij(lon, lat, delta_x, delta_y, NX, NY):

    ''''
    '''
    tlon = int((lon+180.0)/delta_x+1.5)

    tlat = int((lat + 90.0) / delta_y + 1.5)

    if( tlon > NX ):
        tlon = tlon-NX

    return tlon, tlat

def co2_interpolation(gc_pres, obspack_pres, gc_co2_value):
    '''

    :param gc_pres:
    :param obspack_pres:
    :return:
    '''

    if obspack_pres >= gc_pres[0]:
        inter_co2_value = gc_co2_value[0]
    else:
        for verith in range(gc_pres.shape[0]):
            if obspack_pres >= gc_pres[verith]:
                inter_co2_value = gc_co2_value[verith-1]*(gc_pres[verith-1]-obspack_pres)/(gc_pres[verith-1]-gc_pres[verith]) \
                                  + gc_co2_value[verith]*(obspack_pres-gc_pres[verith])/(gc_pres[verith-1]-gc_pres[verith])
                break

    return inter_co2_value

def Barometri_formula(pressure):
    '''

    :param pressure:
    :return:
    '''

    p0 = 101325
    cp = 1004.68506
    g = 9.80665
    M = 0.02896968
    R0 = 8.314462618
    L = 0.00976
    T0 = 288.16

    h = (1-np.exp(np.log(pressure/p0)*R0/(cp*M)))*(cp*T0)/g
    return h


lat_num = 91
lon_num = 144
delta_x = 2.5
delta_y = 2


gc_level = 47
obspack_units = "seconds since 1970-01-01T00:00:00Z"
begin_time_str = '2018-01-01 00:00:00'
end_time_str = '2019-01-01 00:00:00'
end_time_date = datetime.datetime.strptime(end_time_str, '%Y-%m-%d %H:%M:%S')
begin_time_date = datetime.datetime.strptime(begin_time_str, '%Y-%m-%d %H:%M:%S')
end_time_num_sin1970s = date2num(end_time_date, obspack_units)
begin_time_num_sin1970s = date2num(begin_time_date, obspack_units)

#sinal whether have different assimilation data
NADIR = 0
NADIR_PLUS_GLINTLAND = 0
NADIR_PLUS_GLINTLAND_GLINTOCEAN = 0
#reading geos-chem file
gc_file_nums = np.zeros((end_time_date-begin_time_date).days*24*3)
gc_read_time_date = begin_time_date
gc_num = 0
while(gc_read_time_date != end_time_date):
    # gc_ori_dir = '/home/huoxiao/n1s1/huoxiao_share/GEOSChem_Simulation_v2/GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
    #         .format(begin_time_date.year, begin_time_date.month, begin_time_date.day, begin_time_date.hour, begin_time_date.minute)
    # gc_na_plus_glint_dir = '/home/huoxiao/n1s1/huoxiao_share/GEOSChem_Simulation_v2/GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
    #         .format(begin_time_date.year, begin_time_date.month, begin_time_date.day, begin_time_date.hour, begin_time_date.minute)
    # gc_met_dir = '/home/huoxiao/n1s1/huoxiao_share/GEOSChem_Simulation_v2/GEOSChem.StateMet.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
    #     .format(begin_time_date.year, begin_time_date.month, begin_time_date.day, begin_time_date.hour, begin_time_date.minute)
    #
    # gc_ori_nc_data = Dataset(gc_ori_dir, 'r')
    # gc_gc_na_plus_glint_nc_data = Dataset(gc_na_plus_glint_dir, 'r')
    # gc_met_nc_data = Dataset(gc_met_dir, 'r')
    gc_file_nums[gc_num] = date2num(gc_read_time_date, obspack_units)
    gc_read_time_date += datetime.timedelta(minutes=60)
    gc_num += 1


obspack_file_dir = '/share/nas1_share1/huoxiao_share/CO2_data/Observation_data/obspack_co2_1_CARBONTRACKER_CT2019_2020-01-16/data/nc'
pbl_file_dir = '/home/huoxiao/GEOS_FP/GEOS_FP_ADJOINT/data/GEOS_2x2.5/MERRA2'

os.chdir(obspack_file_dir)
obspack_file_name_list = os.listdir()

num = 0
gc_ori_flask_bias_list = []
gc_na_plus_glint_flask_bias_list = []

latitude_list = []
longitude_list = []


obspack_dic_list = []
for obspackith in range(len(obspack_file_name_list)):


    obspack_nc_data = Dataset(obspack_file_name_list[obspackith])

    if(obspack_nc_data.provider_1_affiliation == 'National Oceanic and Atmospheric Administration'):

        if (fnmatch(obspack_file_name_list[obspackith], '*_aircraft-*')):


            if obspack_nc_data.variables['time'][-1] >= end_time_num_sin1970s and obspack_nc_data.variables['time'][0] <= begin_time_num_sin1970s:

                obspack_pressure_list = []
                obspack_co2_list = []
                gc_ori_co2_list = []
                gc_na_plus_glint_co2_list = []
                obspack_time_list = []
                obspack_longitude_list = []
                obspack_latitude_list = []
                lat_index_list = []
                lon_index_list = []
                blh_list = []

                for dataith in range(obspack_nc_data.variables['time'][:].shape[0]):
                    obspack_time_value = obspack_nc_data.variables['time'][dataith]

                    if (obspack_nc_data.variables['time'][dataith] <= end_time_num_sin1970s and
                            obspack_nc_data.variables['time'][dataith] >= begin_time_num_sin1970s and
                            obspack_nc_data.variables['obs_flag'][dataith] == 1):

                        latitude = obspack_nc_data.variables['latitude'][dataith]
                        longitude = obspack_nc_data.variables['longitude'][dataith]
                        obspack_time_value = obspack_nc_data.variables['time'][dataith]
                        obspack_pres_value = obspack_nc_data.variables['pressure'][dataith]


                        gc_file_num_lt_obst_plus_30m = gc_file_nums[np.where(gc_file_nums <= (obspack_time_value + 30 * 60))]
                        gc_file_num = gc_file_num_lt_obst_plus_30m[np.where(gc_file_num_lt_obst_plus_30m >= (obspack_time_value - 30 * 60))]
                        gc_mean_ori_co2_value = np.zeros((gc_level, lat_num, lon_num))
                        gc_mean_na_plus_glint_co2_value = np.zeros((gc_level, lat_num, lon_num))
                        gc_mean_pres_value = np.zeros((gc_level, lat_num, lon_num))
                        pbl_mean_value =  np.zeros((lat_num, lon_num))
                        for gcith in range(gc_file_num.shape[0]):
                            gc_file_date = num2date(gc_file_num[gcith], obspack_units)
                            gc_ori_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                                .format(gc_file_date.year, gc_file_date.month, gc_file_date.day, gc_file_date.hour,
                                        gc_file_date.minute)
                            gc_na_plus_glint_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Assimilation_Litev9_Nadir_Obsm_v2/GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                                .format(gc_file_date.year, gc_file_date.month, gc_file_date.day, gc_file_date.hour,
                                        gc_file_date.minute)
                            gc_met_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/GEOSChem.StateMet.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                                .format(gc_file_date.year, gc_file_date.month, gc_file_date.day, gc_file_date.hour,
                                        gc_file_date.minute)
                            pbl_dir = pbl_file_dir+'/{:0>4d}/{:0>2d}/MERRA2.{:0>4d}{:0>2d}{:0>2d}.A1.2x25.nc4' \
                                .format(gc_file_date.year, gc_file_date.month,
                                        gc_file_date.year, gc_file_date.month, gc_file_date.day)

                            gc_ori_nc_data = Dataset(gc_ori_dir, 'r')

                            gc_na_plus_glint_nc_data = Dataset(gc_na_plus_glint_dir, 'r')
                            gc_met_nc_data = Dataset(gc_met_dir, 'r')

                            pbl_nc_data = Dataset(pbl_dir, 'r')

                            gc_mean_ori_co2_value += gc_ori_nc_data.variables['SpeciesConc_CO2'][0, :, :, :] / gc_file_num.shape[0]
                            gc_mean_na_plus_glint_co2_value += gc_na_plus_glint_nc_data.variables['SpeciesConc_CO2'][0,:, :, :]/gc_file_num.shape[0]
                            gc_mean_pres_value += gc_met_nc_data.variables['Met_PMID'][0, :, :, :] / gc_file_num.shape[0]
                            pbl_mean_value += pbl_nc_data.variables['PBLH'][gc_file_date.hour, :, :]/ gc_file_num.shape[0]

                            gc_ori_nc_data.close()
                            gc_met_nc_data.close()
                            pbl_nc_data.close()
                        # obspack_time_lt_end_time = obspack_nc_data.variables['time'][np.where(obspack_nc_data.variables['time'][:]<end_time_num_sin1970s)]
                        # obspack_time_ht_begin_time = obspack_time_lt_end_time[np.where(obspack_time_lt_end_time>begin_time_num_sin1970s)]
                        # print(num2date(obspack_nc_data.variables['time'][-1], obspack_units))
                        # print(obspack_file_name_list[obspackith])
                        # print(obspack_time_ht_begin_time.shape)
                        # num += 1
                        # latitude_list.append(latitude)
                        # longitude_list.append(longitude)

                        lon_index, lat_index = get_ij(longitude, latitude, delta_x, delta_y, lon_num, lat_num)
                        # fig = plt.figure()
                        # ax = fig.add_subplot()
                        # ax.plot(gc_mean_ori_co2_value[:, lat_index, lon_index]* 10 ** 6, gc_mean_pres_value[:, lat_index, lon_index] * 100)
                        # ax.plot(gc_mean_na_plus_glint_co2_value[:, lat_index, lon_index]* 10 ** 6,gc_mean_pres_value[:, lat_index, lon_index] * 100)
                        # print(sfdf)

                        gc_ori_inter_co2_value = co2_interpolation(gc_mean_pres_value[:, lat_index, lon_index] * 100,obspack_pres_value, \
                                                                   gc_mean_ori_co2_value[:, lat_index, lon_index])
                        gc_na_plus_glint_inter_co2_value = co2_interpolation(gc_mean_pres_value[:, lat_index, lon_index] * 100, obspack_pres_value, \
                            gc_mean_na_plus_glint_co2_value[:, lat_index, lon_index])

                        #append
                        obspack_pressure_list.append(obspack_pres_value)
                        blh_list.append(pbl_mean_value[lat_index, lon_index])
                        obspack_co2_list.append(obspack_nc_data.variables['value'][dataith])
                        gc_ori_co2_list.append(gc_ori_inter_co2_value)
                        gc_na_plus_glint_co2_list.append(gc_na_plus_glint_inter_co2_value)


                        obspack_longitude_list.append(longitude)
                        obspack_latitude_list.append(latitude)
                        obspack_time_list.append(str(num2date(obspack_time_value, obspack_units)))
                        lat_index_list.append(lat_index)
                        lon_index_list.append(lon_index)

                if(len(obspack_time_list)==0):
                    continue
                #for dataith
                print(obspack_file_name_list[obspackith], obspack_nc_data.site_name)
                obspack_dic = {}
                obspack_dic.update({'obspack_time': np.array(obspack_time_list)})
                obspack_dic.update({'obspack_longitude': np.array(obspack_longitude_list)})
                obspack_dic.update({'obspack_latitude': np.array(obspack_latitude_list)})
                obspack_dic.update({'gc_lat_index': np.array(lat_index)})
                obspack_dic.update({'gc_lon_index': np.array(lon_index)})

                obspack_dic.update({'obspack_pressure': np.array(obspack_pressure_list)})
                obspack_dic.update({'obspack_co2': np.array(obspack_co2_list)})
                obspack_dic.update({'gc_ori_co2': np.array(gc_ori_co2_list)})
                obspack_dic.update({'gc_na_plus_glint_co2': np.array(gc_na_plus_glint_co2_list)})

                obspack_dic.update({'obspack_blh': np.array(blh_list)})
                obspack_dic_list.append(obspack_dic)

aircraft_out_dic = {}
pddata = pd.DataFrame(aircraft_out_dic)
for dicith in range(len(obspack_dic_list)):
    pddata = pd.concat([pddata, pd.DataFrame(obspack_dic_list[dicith])], axis=1)
with pd.ExcelWriter('/home/huoxiao/paper_data/gc_aircraft_obs_flag_gc''.xlsx', \
                        mode='w') as writer:
    pddata.to_excel(writer, sheet_name='TCCON SITE DATA', engine='xlsxwriter')


