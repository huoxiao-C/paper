import sys
from netCDF4 import Dataset, num2date, date2num
import datetime
import subprocess
import numpy as np
import calendar
import os
import time
import pandas as pd

def pressure_weigt_function(pres, surfp, level):
    ''''

    '''
    h = np.zeros(level)
    for i in range(level):
        if i == 0:
            h[i] = abs(-pres[i] + (pres[i + 1] - pres[i]) / np.log(pres[i + 1] / pres[i]))/surfp
        elif i == level-1:
            h[i] = abs(pres[i]-(pres[i] - pres[i - 1]) / np.log(pres[i] / pres[i - 1]))/surfp
        else:
            h[i] = abs(-pres[i] + (pres[i + 1] - pres[i]) / np.log(pres[i + 1] / pres[i])+pres[i]-(pres[i] - pres[i - 1]) / np.log(pres[i] / pres[i - 1]))/surfp

    return h

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

def convert_wmf2dmf(wmf, H2o_wmf):
    '''

    :param wmf:
    :param H2o_wmf:
    :return:
    '''

    return wmf/(1-H2o_wmf)

def avk_co2_mean(avk_index, avk_co2):
    '''

    :param avk_index:
    :param avk_co2:
    :return:
    '''
    num = avk_index.shape[0]
    avk_co2_mean = np.zeros(avk_co2.shape[0])
    for i in range(num):
        avk_co2_mean[:] += avk_co2[:, int(avk_index[i])]/num
    return avk_co2_mean


def tccon_concentration_mean(concentration, level):
    ''''
    '''
    mean_concentration = np.zeros(level-1)
    for i in range(level):
        mean_concentration = np.mean(concentration[i:i+1])
    return mean_concentration

def search_avk_index(ak_zenith, asza_deg):
    ''''
    '''
    avk_index = np.zeros(asza_deg.shape[0])

    avk_index[np.where(asza_deg<10)] = 0
    avk_index[np.where(asza_deg > 85)] = 15
    asza_lt10 = np.where(asza_deg < 10)
    asza_ht85 = np.where(asza_deg > 85)


    asza_deg[np.where(asza_deg%5>=2.5)] = np.int32(np.floor(asza_deg[np.where(asza_deg%5>=2.5)]/5)*5+5)
    asza_deg[np.where(asza_deg%5<2.5)] = np.int32(np.floor(asza_deg[np.where(asza_deg % 5 < 2.5)]/5)*5)
    asza_deg[asza_lt10] = -999
    asza_deg[asza_ht85] = -999

    for i in range(ak_zenith.shape[0]):
        avk_index[np.where(ak_zenith[i]==asza_deg)] = i
    return avk_index

def gc_interp_tccon(gc_co2_native, gc_level, gc_pres, tccon_level, tccon_pres):
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

def convert_co2_to_xco2(gc_co2, tccon_avk_co2, tccon_prior_co2, tccon_pres, tccon_level, \
                        tccon_fH2o, constant_g, Mh2o, Mair):
    '''

    :param gc_co2:
    :param tccon_avk_co2:
    :param tccon_prior_co2:
    :param tccon_pres:
    :return:
    '''

    tmp_gc_co2 = np.zeros(tccon_level-1)
    tmp_prior_co2 = np.zeros(tccon_level-1)
    tmp_delta_pres = np.zeros(tccon_level-1)
    tmp_pres =  np.zeros(tccon_level-1)
    tmp_avk_co2 = np.zeros(tccon_level-1)
    tmp_fH2o = np.zeros(tccon_level-1)
    tmp_constant_g = np.zeros(tccon_level-1)
    for level in range(tccon_level-1):
        tmp_gc_co2[level] = np.mean(gc_co2[level:level+2])
        tmp_prior_co2[level] = np.mean(tccon_prior_co2[level:level + 2])
        tmp_delta_pres[level] = tccon_pres[level+1]-tccon_pres[level]
        tmp_pres[level] = np.mean(tccon_pres[level:level+2])
        tmp_avk_co2[level] = np.mean(tccon_avk_co2[level:level+2])
        tmp_fH2o[level] = np.mean(tccon_fH2o[level:level+2])
        tmp_constant_g = np.mean(constant_g[level:level+2])

    # VC_gc_co2 = np.sum(tmp_gc_co2*tmp_avk_co2*tmp_delta_pres/(tmp_constant_g*Mair*(1+tmp_fH2o*(Mh2o/Mair))))
    # VC_prior_co2 = np.sum(tmp_prior_co2*tmp_avk_co2*tmp_delta_pres/(tmp_constant_g*Mair*(1+tmp_fH2o*(Mh2o/Mair))))
    # VC_air = np.sum(tmp_delta_pres/(tmp_constant_g * Mair * (1 + tmp_fH2o * (Mh2o / Mair))))

    # tccon_prior_xco2 = np.sum(tmp_prior_co2*tmp_delta_pres/(tmp_constant_g*Mair*(1+tmp_fH2o*(Mh2o/Mair))))/VC_air

    # gc_xco2  = tccon_prior_xco2+(VC_gc_co2-VC_prior_co2)/VC_air

    #pressure weighting function
    h = pressure_weigt_function(tmp_pres, tmp_pres[0], tccon_level-1)
    tccon_prior_xco2 = np.sum(h*tmp_prior_co2)

    gc_xco2 = tccon_prior_xco2+np.sum(h*(tmp_gc_co2-tmp_prior_co2))
    return gc_xco2

def RMSE(pred, obs):
    '''

    :param pred:
    :param obs:
    :return:
    '''
    return np.sqrt(((pred-obs)**2).mean())

pbtime = time.time()

tccon_dir_path = '/share/nas1_share1/huoxiao_share/CO2_data/Observation_data/tccon/tccon_all/'
gcrspc_pri_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/'
gcmet_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/'
gcrspc_pos_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Assimilation_Litev9_Nadir_Obsm_v2/'
#tccon_omit_all_name = ('ae', 'an', 'br', 'bu', 'ci', 'db', 'df', 'et', 'eu', 'gm', 'iz', 'js', \
#                    'ka', 'll', 'oc', 'or', 'pa', 'pr', 'ra', 'rj', 'so', 'sp', 'tk', 'wg')
#tccon_all_name = ['Ascension','Anmyeondo','Breman', 'Burgos', 'Caltech', 'Darwin', 'Dryden', 'East', 'Eureka', 'Garmisch', \
                  #'Izana', 'Saga', 'Kalsruche', 'Lauder','Lamont','Orleans', 'Park Falls', 'Paris', 'Reunion', 'Rikubetsu', \
                  #'Sodankula', 'Ny Alesund', 'Tsukuba', 'Wollongong']

tccon_all_name = ['Ascension Island', 'Anmyeondo, Korea', 'Sites/Bialystok', 'Bremen,Germany', 'Burgos, Philippines', 'Caltech, USA', \
                  'Darwin,Australia', 'Dryden, USA', 'East Trout Lake, Canada', 'Eureka,Canada', 'Four Corners, USA', 'Garmisch,Germany', \
                  'Hefei, China', 'Indianapolis, IN, USA', 'Izana,Tenerife', 'JPL, Pasadena, CA, USA', 'JPL, Pasadena, CA, USA', 'Saga, Japan',\
                  'Karlsruhe, Germany', 'Lauder,New Zealand', 'Lauder,New Zealand', 'Lauder,New Zealand', 'Manaus, Brazil', 'Lamont,OK (USA)', \
                  'Orleans,France', 'Park Falls,WI (USA)', 'Paris, France', 'Reunion Island', 'Rikubetsu, Japan ', 'Sodankyla, Finland ', \
                  'Ny Alesund,Spitsbergen', 'Tsukuba,Japan', 'Wollongong,Australia', 'Zugspitze,Germany']


temp = subprocess.check_output('ls /share/nas1_share1/huoxiao_share/CO2_data/Observation_data/tccon/tccon_all/', shell=True)
tccon_time_units = 'days since 1970-01-01 00:00:00'
#convert time to number days since 1970-01-01 00:00:00
co2_begin_time_str = "2018-01-01 00:00:00"
co2_end_time_str = "2019-01-01 00:00:00"

#
#list tccon file
command_output = subprocess.check_output('ls '+tccon_dir_path, shell=True)
tccon_filenames = files_string2list(command_output.decode())




#geos-chem data read

#parameter
Mair = 29.0 * 10 ** -3 / (6.022 * 10 ** 23)
Mh2o = 18.0 ** -3 / (6.022 * 10 ** 23)
gc_level = 47
gc_lat = 91
gc_lon = 144
tccon_level = 71
g = 9.8
# gc_file_nums = np.zeros(((co2_end_time_date-co2_begin_time_date).days*24, gc_level, gc_lat, gc_lon))




#oco2
# satellite
# oco2_begin_time = co2_begin_time_date
# oco2_end_time = co2_end_time_date
# oco2_all_xco2_list = []
# oco2_all_lon_list = []
# oco2_all_lat_list = []
#
# while (oco2_begin_time!=oco2_end_time):
#     tmp_xco2_file = np.array([0])
#     tmp_lon_file = np.array([0])
#     tmp_lat_file = np.array([0])
#
#     oco2_day_end_time = oco2_begin_time+datetime.timedelta(days=1)
#     while (oco2_begin_time != oco2_day_end_time):
#         filename = '/home/huoxiao/n1s1/huoxiao_share/CO2_data/Observation_data/satellite_on_gc_2x25grid/OCO2_CO2_Litev9_ND_plus_GlintLand/oco2_' + oco2_begin_time.strftime(
#             '%Y%m%d%H%M%S') + '.nc'
#         if not subprocess.call('test -f ' + filename, shell=True):
#             print(filename)
#             oco2_nc_data = Dataset(filename, 'r')
#             tmp_xco2_file = np.concatenate((tmp_xco2_file, oco2_nc_data.variables['xco2'][:]), axis=0)
#             tmp_lat_file = np.concatenate((tmp_lat_file, oco2_nc_data.variables['latitude'][:]), axis=0)
#             tmp_lon_file = np.concatenate((tmp_lon_file, oco2_nc_data.variables['longitude'][:]), axis=0)
#         oco2_begin_time = oco2_begin_time + datetime.timedelta(minutes=20)
#     if(tmp_xco2_file[1::].shape[0]==0):
#         oco2_all_xco2_list.append(np.array([-1]))
#         oco2_all_lon_list.append(np.array([-1]))
#         oco2_all_lat_list.append(np.array([-1]))
#
#     oco2_all_xco2_list.append(tmp_xco2_file[1::])
#     oco2_all_lon_list.append(tmp_lon_file[1::])
#     oco2_all_lat_list.append(tmp_lat_file[1::])





#sinal whether have different assimilation data
NADIR = 1
NADIR_PLUS_GLINTLAND = 0
NADIR_PLUS_GLINTLAND_GLINTOCEAN = 0
oco2_obs_flag_array = np.array([NADIR, NADIR_PLUS_GLINTLAND, NADIR_PLUS_GLINTLAND_GLINTOCEAN])
oco2_obs_flag_list = ['nd', 'nd_plus_glintland', 'nd_plus_glintland_plus_glintocean']
for oco2_obsith in range(oco2_obs_flag_array.shape[0]):
    bias_array_all = np.zeros((6, 50))
    rmse_array = np.zeros((6, 50))
    all_pri_data_hour = np.zeros((3, 50, 100000))
    all_pos_data_hour = np.zeros((3, 50, 100000))
    all_tccon_data_hour = np.zeros((3, 50, 100000))
    all_time_data_hour = np.empty(shape=(3, 50, 100000), dtype='<U30')
    bias_array_hour_num_array = np.zeros((3, 50))

    co2_begin_time_date = datetime.datetime.strptime(co2_begin_time_str, '%Y-%m-%d %H:%M:%S')
    co2_end_time_date = datetime.datetime.strptime(co2_end_time_str, '%Y-%m-%d %H:%M:%S')

    tccon_filename_list = []
    site_list = []
    if oco2_obs_flag_array[oco2_obsith] == 0:
        continue


    for tccon_i in range(len(tccon_filenames)):


        tccon_path = tccon_dir_path+tccon_filenames[tccon_i]
        print(tccon_i, tccon_path)
        tccon_nc_data = Dataset(tccon_path, 'r')

        gc_pos_xco2_array_hour = np.zeros(100000)
        gc_pri_xco2_array_hour = np.zeros(100000)
        tccon_xco2_array_hour = np.zeros(100000)
        obs_time_array_hour = np.empty(shape=(100000), dtype='<U30')

        bias_array_hour_num = 0
        tccon_dic_list = []



        co2_begin_time_date = datetime.datetime.strptime(co2_begin_time_str, '%Y-%m-%d %H:%M:%S')
        co2_end_time_date = datetime.datetime.strptime(co2_end_time_str, '%Y-%m-%d %H:%M:%S')
        while(co2_begin_time_date!=co2_end_time_date):

            tccon_date_num = tccon_nc_data.variables['time'][:].shape[0]
            tccon_second_num = tccon_nc_data.variables['prior_date'][:].shape[0]

            co2_begin_time_num = date2num(co2_begin_time_date-datetime.timedelta(minutes=30), tccon_time_units)
            co2_end_time_num = date2num(co2_begin_time_date+datetime.timedelta(minutes=30), tccon_time_units)

            #time
            tccon_time_lt_end_time = tccon_nc_data.variables['time'][np.where(tccon_nc_data.variables['time'][:]<=co2_end_time_num)]
            tccon_time_ht_end_time = tccon_time_lt_end_time[np.where(tccon_time_lt_end_time>co2_begin_time_num)]

            if tccon_time_ht_end_time.shape[0] == 0:
                co2_begin_time_date += datetime.timedelta(hours=1)

                continue

            tccon_xco2_lt_end_time = tccon_nc_data.variables['xco2_ppm'][np.where(tccon_nc_data.variables['time'][:] <= co2_end_time_num)]
            tccon_asza_deg_lt_end_time = tccon_nc_data.variables['asza_deg'][np.where(tccon_nc_data.variables['time'][:] <= co2_end_time_num)]

            tccon_xco2_ht_end_time = tccon_xco2_lt_end_time[np.where(tccon_time_lt_end_time>=co2_begin_time_num)]
            tccon_asza_deg_ht_end_time = tccon_asza_deg_lt_end_time[np.where(tccon_time_lt_end_time>=co2_begin_time_num)]


            #search avk 2th dimension index
            avk_index = search_avk_index(tccon_nc_data.variables['ak_zenith'][:], tccon_asza_deg_ht_end_time)
            tccon_avk_co2 =  avk_co2_mean(avk_index, tccon_nc_data.variables['ak_co2'][:, :])

            #prior date
            #floor time num
            co2_begin_time_num = np.round(co2_begin_time_num)
            co2_end_time_num = np.round(co2_end_time_num)

            tccon_pdate_lt_end_time = tccon_nc_data.variables['prior_date'][np.where(tccon_nc_data.variables['prior_date'][:]<=co2_end_time_num)]
            tccon_prior_gravity_lt_end_time = tccon_nc_data.variables['prior_gravity'][np.where(tccon_nc_data.variables['prior_date'][:] <= co2_end_time_num)[0], :]
            tccon_prior_Pressure_lt_end_time = tccon_nc_data.variables['prior_Pressure'][np.where(tccon_nc_data.variables['prior_date'][:] <= co2_end_time_num)[0], :]
            tccon_prior_h2o_lt_end_time = tccon_nc_data.variables['prior_h2o'][np.where(tccon_nc_data.variables['prior_date'][:] <= co2_end_time_num)[0], :]
            tccon_prior_co2_lt_end_time = tccon_nc_data.variables['prior_co2'][np.where(tccon_nc_data.variables['prior_date'][:] <= co2_end_time_num)[0], :]

            tccon_pdate_ht_end_time = tccon_pdate_lt_end_time[np.where(tccon_pdate_lt_end_time>=co2_begin_time_num)]
            tccon_prior_gravity_ht_end_time = tccon_prior_gravity_lt_end_time[np.where(tccon_pdate_lt_end_time>=co2_begin_time_num)[0], :]
            tccon_prior_Pressure_ht_end_time = tccon_prior_Pressure_lt_end_time[np.where(tccon_pdate_lt_end_time>=co2_begin_time_num)[0], :]
            tccon_prior_h2o_ht_end_time = tccon_prior_h2o_lt_end_time[np.where(tccon_pdate_lt_end_time>=co2_begin_time_num)[0], :]
            tccon_prior_co2_ht_end_time = tccon_prior_co2_lt_end_time[np.where(tccon_pdate_lt_end_time>=co2_begin_time_num)[0], :]

            #convert to dry air
            tccon_h2o_dry = tccon_prior_h2o_ht_end_time/(1-tccon_prior_h2o_ht_end_time)
            tccon_pri_co2_dry = tccon_prior_co2_ht_end_time/(1-tccon_prior_h2o_ht_end_time)


            if tccon_nc_data.variables['prior_date'][np.where(tccon_nc_data.variables['prior_date'][:]==co2_begin_time_num)].shape[0] == 0:
                co2_begin_time_date += datetime.timedelta(hours=1)

                continue

            tccon_dic = {}

            tccon_dic.update({'xco2': tccon_xco2_ht_end_time*10**-6})
            tccon_dic.update({'prior_co2': tccon_pri_co2_dry*10**-6})
            tccon_dic.update({'avk_co2': tccon_avk_co2})
            tccon_dic.update({'gravity': tccon_prior_gravity_ht_end_time})
            tccon_dic.update({'prior_h2o': tccon_h2o_dry})
            tccon_dic.update({'prior_pressure': tccon_prior_Pressure_ht_end_time})
            tccon_dic.update({'lon': tccon_nc_data.variables['long_deg'][0]})
            tccon_dic.update({'lat': tccon_nc_data.variables['lat_deg'][0]})


            tccon_filename_list.append(tccon_filenames[tccon_i][0:2])

            tccon_dic_list.append(tccon_dic)


            # GEOS-CHEM file read
            # gc



            # search lon lat index

            # set lat and lon span
            if tccon_nc_data.variables['lat_deg'][0] > -25:
                lat_span = 2.5
                lon_span = 5
            else:
                lat_span = 5
                lon_span = 60
            # calendar.monthrange(2016, month)[1]+1
            #read gc file
            gcrspc_pri_filename = 'GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4' \
                .format(co2_begin_time_date.year, co2_begin_time_date.month, co2_begin_time_date.day, co2_begin_time_date.hour)
            gcmet_filename = 'GEOSChem.StateMet.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4' \
                .format(co2_begin_time_date.year, co2_begin_time_date.month, co2_begin_time_date.day, co2_begin_time_date.hour)
            gcrspc_pos_filename = 'GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}00z.nc4' \
                    .format(co2_begin_time_date.year, co2_begin_time_date.month, co2_begin_time_date.day, co2_begin_time_date.hour)
            gcrspc_pos_path = gcrspc_pos_dir+gcrspc_pos_filename
            gcrspc_pri_path = gcrspc_pri_dir+gcrspc_pri_filename
            gcmet_path = gcmet_dir+gcmet_filename

            gcspc_pri_nc_data = Dataset(gcrspc_pri_path, 'r')
            gcmet_nc_data = Dataset(gcmet_path, 'r')
            gcspc_pos_nc_data = Dataset(gcrspc_pos_path, 'r')

            #search lon lat index
            gc_lon_lt_tccon_lon = gcmet_nc_data.variables['lon'][np.where(gcmet_nc_data.variables['lon'][:] <= (tccon_dic['lon'] + lon_span))]
            gc_lon_ht_tccon_lon = gc_lon_lt_tccon_lon[gc_lon_lt_tccon_lon >= (tccon_dic['lon'] - lon_span)]

            gc_lat_lt_tccon_lat = gcmet_nc_data.variables['lat'][np.where(gcmet_nc_data.variables['lat'][:] <= (tccon_dic['lat'] + lat_span))]
            gc_lat_ht_tccon_lat = gc_lat_lt_tccon_lat[gc_lat_lt_tccon_lat >= (tccon_dic['lat'] - lat_span)]


            gc_pos_co2 = gcspc_pos_nc_data.variables['SpeciesConc_CO2'][0, :, :, :]

            gc_pri_co2 = gcspc_pri_nc_data.variables['SpeciesConc_CO2'][0, :, :, :]
            gc_met_pres = gcmet_nc_data.variables['Met_PMID'][0, :, :, :]

            #prior
            gc_lon_num = gc_lon_ht_tccon_lon.shape[0]
            gc_lat_num = gc_lat_ht_tccon_lat.shape[0]
            gc_co2_native_ori = np.zeros(gc_level)
            gc_co2_native_ass = np.zeros(gc_level)
            gc_pres = np.zeros(gc_level)
            for i in range(gc_lon_num):
                for j in range(gc_lat_num):
                    gc_co2_native_ori += (gc_pri_co2[:, np.where(gc_lat_ht_tccon_lat[j] == gcmet_nc_data.variables['lat'][:])[0], \
                                          np.where(gc_lon_ht_tccon_lon[i] == gcmet_nc_data.variables['lon'][:])[0]] \
                                          / (gc_lon_num * gc_lat_num)).reshape(gc_level)
                    gc_pres += (gc_met_pres[:, np.where(gc_lat_ht_tccon_lat[j] == gcmet_nc_data.variables['lat'][:])[0], \
                                np.where(gc_lon_ht_tccon_lon[i] == gcmet_nc_data.variables['lon'][:])[0]] \
                                / (gc_lon_num * gc_lat_num)).reshape(gc_level)
                    gc_co2_native_ass += (gc_pos_co2[:, np.where(gc_lat_ht_tccon_lat[j] == gcmet_nc_data.variables['lat'][:])[0],
                                np.where(gc_lon_ht_tccon_lon[i] == gcmet_nc_data.variables['lon'][:])[0]] \
                                / (gc_lon_num * gc_lat_num)).reshape(gc_level)

            gc_co2_ori = gc_interp_tccon(gc_co2_native_ori, gc_level, gc_pres, tccon_level, np.mean(np.array(tccon_dic['prior_pressure']), axis=0))
            gc_co2_ass = gc_interp_tccon(gc_co2_native_ass, gc_level, gc_pres, tccon_level,np.mean(np.array(tccon_dic['prior_pressure']), axis=0))


            gc_xco2_ori = convert_co2_to_xco2(gc_co2_ori, tccon_dic['avk_co2'], \
                                          np.mean(tccon_dic['prior_co2'], axis=0),       \
                                          np.mean(tccon_dic['prior_pressure'], axis=0), \
                                          tccon_level, np.mean(tccon_dic['prior_h2o'], axis=0), \
                                          np.mean(tccon_dic['gravity'], axis=0), Mh2o, Mair)

            gc_xco2_ass = convert_co2_to_xco2(gc_co2_ass, tccon_dic['avk_co2'], \
                                          np.mean(tccon_dic['prior_co2'], axis=0),       \
                                          np.mean(tccon_dic['prior_pressure'], axis=0), \
                                          tccon_level, np.mean(tccon_dic['prior_h2o'], axis=0), \
                                          np.mean(tccon_dic['gravity'], axis=0), Mh2o, Mair)


            #increment

            gc_pri_xco2_array_hour[bias_array_hour_num] = gc_xco2_ori
            gc_pos_xco2_array_hour[bias_array_hour_num] = gc_xco2_ass
            tccon_xco2_array_hour[bias_array_hour_num] = np.mean(np.array(tccon_dic['xco2']))
            obs_time_array_hour[bias_array_hour_num] = datetime.datetime.strftime(co2_begin_time_date, '%Y-%m-%d %H:%M:%S')

            bias_array_hour_num += 1
            co2_begin_time_date += datetime.timedelta(hours=1)
            gcmet_nc_data.close()
            gcspc_pri_nc_data.close()
            gcspc_pos_nc_data.close()
            # print('ori', gc_xco2_ori-np.mean(np.array(tccon_dic['xco2'])), gc_xco2_ori, np.mean(np.array(tccon_dic['xco2'])))
            # print('ass', gc_xco2_ass-np.mean(np.array(tccon_dic['xco2'])), gc_xco2_ass, np.mean(np.array(tccon_dic['xco2'])))

        #judge
        if(bias_array_hour_num ==0):
            continue
        #bias
        site_list.append(tccon_i)


        all_pri_data_hour[oco2_obsith, len(site_list)-1, 0:bias_array_hour_num] = gc_pri_xco2_array_hour[0:bias_array_hour_num]
        all_pos_data_hour[oco2_obsith, len(site_list)-1, 0:bias_array_hour_num] = gc_pos_xco2_array_hour[0:bias_array_hour_num]
        all_tccon_data_hour[oco2_obsith, len(site_list)-1, 0:bias_array_hour_num] = tccon_xco2_array_hour[0:bias_array_hour_num]
        print(gc_pri_xco2_array_hour[0:bias_array_hour_num], 'number', bias_array_hour_num)
        all_time_data_hour[oco2_obsith, len(site_list)-1, 0:bias_array_hour_num] = obs_time_array_hour[0:bias_array_hour_num]
        bias_array_hour_num_array[oco2_obsith, len(site_list)-1] = bias_array_hour_num
        tccon_nc_data.close()

    tccon_cur_name = []
    tccon_out_dic = {}
    pddata = pd.DataFrame(tccon_out_dic)
    for i in range(len(site_list)):
        site_name = tccon_all_name[site_list[i]]
        obs_num = int(bias_array_hour_num_array[oco2_obsith, i])
        tccon_out_dic.update({site_name + ' time': all_time_data_hour[oco2_obsith, i, 0:obs_num]})
        tccon_out_dic.update({site_name+' prior': all_pri_data_hour[oco2_obsith, i, 0:obs_num]})
        tccon_out_dic.update({site_name + ' posterior': all_pos_data_hour[oco2_obsith, i, 0:obs_num]})
        tccon_out_dic.update({site_name + ' tccon': all_tccon_data_hour[oco2_obsith, i, 0:obs_num]})
        pddata = pd.concat([pddata, pd.DataFrame(tccon_out_dic)], axis=1)

        tccon_out_dic = {}

    with pd.ExcelWriter('/home/huoxiao/paper_data/tccon_gc_single_site_'+oco2_obs_flag_list[oco2_obsith]+'.xlsx', \
                        mode='w') as writer:

        pddata.to_excel(writer, sheet_name='TCCON SITE DATA', engine='xlsxwriter')

petime = time.time()
print('program running time:', petime-pbtime)
