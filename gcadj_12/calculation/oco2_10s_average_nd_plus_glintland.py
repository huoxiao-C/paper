from netCDF4 import Dataset, date2num, num2date
import os, time
import subprocess
import ctypes
from matplotlib import pyplot as plt
import numpy as np
import datetime
import warnings
from decimal import Decimal


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


def oco2_weighting_calculation(xco2_day,  xco2_uncert_day, lat_day, lon_day, co2_prior_day,  xco2_prior_day, xco2_avk_day, pres_day, \
            surfp_day, pressure_weight_day, time_day, start_index, last_index, is_10s):
    '''
    浜屾哀鍖栫⒊寤撶嚎娴撳害, 鏂瑰樊1s璁＄畻
    :param oco2:
    :param first_index:
    :param last_index:
    :return:
    '''



    N = last_index-start_index

    #delte nan data
    for i in range(start_index, last_index):
        if(np.isnan(xco2_uncert_day[i])):
            print('nan true')

            if i>start_index:
                xco2_uncert_day[i] = xco2_uncert_day[i-1]
            else:
                xco2_uncert_day[i] = xco2_uncert_day[i+1]

    xco2_inverse_sigma_square = 1/np.square(xco2_uncert_day[start_index:last_index])

    xco2bar = np.sum(xco2_day[start_index:last_index]*xco2_inverse_sigma_square/np.sum(xco2_inverse_sigma_square))
    latbar = np.sum(lat_day[start_index:last_index]*xco2_inverse_sigma_square/np.sum(xco2_inverse_sigma_square))
    co2_priorbar = np.dot(co2_prior_day[start_index:last_index].transpose(), xco2_inverse_sigma_square.transpose())/np.sum(xco2_inverse_sigma_square)
    xco2_priorbar = np.sum(xco2_prior_day[start_index:last_index] * xco2_inverse_sigma_square) / np.sum(xco2_inverse_sigma_square)
    xco2_avkbar = np.dot(xco2_avk_day[start_index:last_index].transpose(), xco2_inverse_sigma_square.transpose()) / np.sum(xco2_inverse_sigma_square)
    presbar = np.dot(pres_day[start_index:last_index].transpose(), xco2_inverse_sigma_square.transpose()) / np.sum(xco2_inverse_sigma_square)
    surfpbar = np.sum(surfp_day[start_index:last_index] * xco2_inverse_sigma_square) / np.sum(xco2_inverse_sigma_square)
    pressure_weightbar = np.dot(pressure_weight_day[start_index:last_index].transpose(), xco2_inverse_sigma_square.transpose()) / np.sum(xco2_inverse_sigma_square)
    timebar = np.sum((time_day[start_index:last_index]-time_day[start_index]) * xco2_inverse_sigma_square)\
              / np.sum(xco2_inverse_sigma_square)+time_day[start_index]

    # #debug
    # debug_unit = 'seconds since 1970-01-01 00:00:00'
    # if timebar-time_day[start_index]<-1:
    #     print(time_day[last_index-1]-time_day[start_index], num2date(time_day[start_index], debug_unit), num2date(time_day[last_index-1], debug_unit))
    # #print(num2date(timebar, debug_unit), num2date(time_day[start_index], debug_unit), num2date(time_day[last_index-1], debug_unit), timebar-time_day[start_index])
    # #end debug


    #base longitude 178 deg
    if(np.max(lon_day[start_index:last_index])>178):
        tmp_lon = np.zeros((last_index-start_index))
        for cur_index in range(start_index, last_index):
            if lon_day[cur_index] < -178:
                tmp_lon[cur_index-start_index] = lon_day[cur_index]+180+2
            else:
                tmp_lon[cur_index - start_index] = lon_day[cur_index] - 178
        tmp_lonbar = np.sum(tmp_lon* xco2_inverse_sigma_square)/np.sum(xco2_inverse_sigma_square)
        if(tmp_lonbar>2):
            lonbar = -180+(tmp_lonbar-2)
        else:
            lonbar = 178+tmp_lonbar

    #base longitude  -178 deg
    elif(np.min(lon_day[start_index])<-178):
        tmp_lon = np.zeros((last_index - start_index))
        for cur_index in range(start_index, last_index):
            if lon_day[cur_index] > 178:
                tmp_lon[cur_index-start_index] = lon_day[cur_index]-180-2
            else:
                tmp_lon[cur_index-start_index] = -180-lon_day[cur_index]

        tmp_lonbar = np.sum(tmp_lon * xco2_inverse_sigma_square) / np.sum(xco2_inverse_sigma_square)

        if(tmp_lonbar<-2):
            lonbar = 180-(-2-tmp_lonbar)
        else:
            lonbar = -178+tmp_lonbar
    else:
        lonbar = np.sum(lon_day[start_index:last_index]*xco2_inverse_sigma_square) / np.sum(xco2_inverse_sigma_square)


    xco2_sigma_IVE = 1/np.sqrt(np.sum(xco2_inverse_sigma_square)/N)
    xco2_floor = np.min(xco2_uncert_day[start_index:last_index])
    xco2_sigma_spread = 0
    if(N!=1):
        xco2_sigma_spread = np.std(xco2_day[start_index:last_index] , ddof=1)

    #variance judgement
    if is_10s:
        xco2_uncertbar = xco2_sigma_IVE
    else:
        xco2_uncertbar = max(xco2_sigma_IVE, xco2_floor, xco2_sigma_spread)


    #reset array
    co2_priorbar = co2_priorbar.transpose()
    xco2_avkbar = xco2_avkbar.transpose()
    presbar = presbar.transpose()
    pressure_weightbar = pressure_weightbar.transpose()



    return xco2bar, xco2_uncertbar, latbar, lonbar, co2_priorbar, xco2_priorbar, xco2_avkbar, presbar, \
           surfpbar, pressure_weightbar, timebar


def oco2_10s_all(xco2_1s, xco2_uncert_1s, lat_1s, lon_1s, co2_prior_1s, xco2_prior_1s,xco2_avk_1s, pres_1s, \
                 surfp_1s, pressure_weight_1s, time_1s, length):
    num = 0
    start_index = 0
    tmp_10s_xco2 = np.zeros(length)
    tmp_10s_xco2_var = np.zeros(length)
    tmp_10s_lat = np.zeros(length)
    tmp_10s_lon = np.zeros(length)
    tmp_10s_co2_apri_profile = np.zeros((length, 20))
    tmp_10s_xco2_apri = np.zeros(length)
    tmp_10s_xco2_avk = np.zeros((length, 20))
    tmp_10s_pres = np.zeros((length, 20))
    tmp_10s_surf = np.zeros(length)
    tmp_10s_xco2_pwf = np.zeros((length, 20))
    tmp_10s_time = np.zeros(length)



    while (True):
        begin_time_num = time_1s[start_index]

        # find time interval index in 10s
        I_flag = True #i exit normally
        for i in range(start_index, length):
            time_num = time_1s[i]
            if (time_num - begin_time_num) < 10:
                continue
            else:
                if i== length-1:
                    I_flag = False
                break

        last_index = i




        if i == length - 1 and I_flag:
            last_index = i+1


        tmp_10s_xco2[num], tmp_10s_xco2_var[num], tmp_10s_lat[num], tmp_10s_lon[num], tmp_10s_co2_apri_profile[num, :], tmp_10s_xco2_apri[num], \
        tmp_10s_xco2_avk[num, :], tmp_10s_pres[num, :], tmp_10s_surf[num], tmp_10s_xco2_pwf[num, :], tmp_10s_time[num] = oco2_weighting_calculation(xco2_1s, \
                xco2_uncert_1s, lat_1s, lon_1s, co2_prior_1s, xco2_prior_1s,xco2_avk_1s, pres_1s, surfp_1s, pressure_weight_1s, \
                time_1s, start_index, last_index, True)

        num += 1
        if i == length - 1 and I_flag:
            break
        # #debug
        # if tmp_10s_lon[num-1]>-50 and tmp_10s_lon[num-1]<-30 and tmp_10s_lat[num-1]>10:
        #     print(time_1s[start_index:last_index], time_1s[last_index-1]-time_1s[start_index], )
        #     print('lon, lat', lon_1s[start_index:last_index], lat_1s[start_index:last_index])
        # #end debug
        start_index = i


    return tmp_10s_xco2, tmp_10s_xco2_var, tmp_10s_lat, tmp_10s_lon, tmp_10s_co2_apri_profile, tmp_10s_xco2_apri, tmp_10s_xco2_avk, tmp_10s_pres, \
           tmp_10s_surf, tmp_10s_xco2_pwf, tmp_10s_time, num

in_oco2_dir = '/share/nas1_share2/zhangqw_share/observation_data/OCO2_DATA/OCO2_L2_Lite_FP.9r/'
out_oco2_dir = '/share/nas1_share1/huoxiao_share/CO2_data/Observation_data/satellite_10s_average/OCO2_Litev9_ND_plus_Glintland/'

oco2_begin_time_str = "2018-12-31 00:00:00"
oco2_end_time_str = "2019-01-01 00:00:00"
oco2_begin_time_date = datetime.datetime.strptime(oco2_begin_time_str, '%Y-%m-%d %H:%M:%S')
oco2_end_time_date = datetime.datetime.strptime(oco2_end_time_str, '%Y-%m-%d %H:%M:%S')

oco2_cur_time_date = oco2_begin_time_date


pro_begin_time = time.time()


#debug
huoxiao_debug_num = 0
warnings.filterwarnings('error')
while(oco2_cur_time_date!=oco2_end_time_date):


    #Determine if the file exists,
    if subprocess.call('test -f ' + in_oco2_dir + '{:0>4d}/oco2_LtCO2_{:0>2d}{:0>2d}{:0>2d}'.format(oco2_cur_time_date.year,oco2_cur_time_date.year%2000, \
        oco2_cur_time_date.month,oco2_cur_time_date.day) + '_B9003r*',shell=True):
        oco2_cur_time_date += datetime.timedelta(days=1)
        continue

    os.chdir(in_oco2_dir)
    command_output = subprocess.check_output("ls  "+ '{:0>4d}/oco2_LtCO2_{:0>2d}{:0>2d}{:0>2d}'.format(oco2_cur_time_date.year,oco2_cur_time_date.year%2000, \
        oco2_cur_time_date.month,oco2_cur_time_date.day) + '_B9003r*',shell=True)

    oco2_filename = files_string2list(command_output.decode())[0]
    oco2_file_path = in_oco2_dir+oco2_filename



    print(oco2_file_path+' is reading')

    in_oco2_nc_data = Dataset(oco2_file_path, 'r')


    # #file is empty
    # if in_oco2.groups['RetrievalResults']['outcome_flag'][:].shape[0] == 0:
    #     continue


    #operation mode: Nadir
    tmp_xco2_day_omeq0 = np.array(in_oco2_nc_data.variables['xco2'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:]==0)]
    tmp_lon_day_omeq0 = np.array(in_oco2_nc_data.variables['longitude'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0)]
    tmp_lat_day_omeq0 = np.array(in_oco2_nc_data.variables['latitude'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0)]
    tmp_co2_prior_day_omeq0 = np.array(in_oco2_nc_data.variables['co2_profile_apriori'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0), :][0]
    tmp_xco2_prior_day_omeq0 = np.array(in_oco2_nc_data.variables['xco2_apriori'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0)]
    tmp_pres_day_omeq0 = np.array(in_oco2_nc_data.variables['pressure_levels'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0), :][0]
    tmp_surfp_day_omeq0 = np.array(in_oco2_nc_data.groups['Retrieval']['psurf'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0)]
    tmp_xco2_uncert_day_omeq0 = np.array(in_oco2_nc_data.variables['xco2_uncertainty'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0)]
    tmp_pressure_weight_day_omeq0 = np.array(in_oco2_nc_data.variables['pressure_weight'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0), :][0]
    tmp_xco2_avk_uncert_day_omeq0 = np.array(in_oco2_nc_data.variables['xco2_averaging_kernel'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0), :][0]
    tmp_time_day_omeq0 = np.array(in_oco2_nc_data.variables['time'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0)]
    tmp_flag_day_omeq0 = np.array(in_oco2_nc_data.variables['xco2_quality_flag'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 0)]



    tmp_xco2_day_omeq0_gd = tmp_xco2_day_omeq0[np.where(tmp_flag_day_omeq0==0)]
    tmp_lon_day_omeq0_gd = tmp_lon_day_omeq0[np.where(tmp_flag_day_omeq0 == 0)]
    tmp_lat_day_omeq0_gd = tmp_lat_day_omeq0[np.where(tmp_flag_day_omeq0 == 0)]
    tmp_co2_prior_day_omeq0_gd = tmp_co2_prior_day_omeq0[np.where(tmp_flag_day_omeq0 == 0), :][0]
    tmp_xco2_prior_day_omeq0_gd = tmp_xco2_prior_day_omeq0[np.where(tmp_flag_day_omeq0 == 0)]
    tmp_pres_day_omeq0_gd = tmp_pres_day_omeq0[np.where(tmp_flag_day_omeq0 == 0), :][0]
    tmp_surfp_day_omeq0_gd = tmp_surfp_day_omeq0[np.where(tmp_flag_day_omeq0 == 0)]
    tmp_xco2_uncert_day_omeq0_gd = tmp_xco2_uncert_day_omeq0[np.where(tmp_flag_day_omeq0 == 0)]
    tmp_pressure_weight_day_omeq0_gd = tmp_pressure_weight_day_omeq0[np.where(tmp_flag_day_omeq0 == 0), :][0]
    tmp_xco2_avk_uncert_day_omeq0_gd = tmp_xco2_avk_uncert_day_omeq0[np.where(tmp_flag_day_omeq0 == 0), :][0]
    tmp_time_day_omeq0_gd = tmp_time_day_omeq0[np.where(tmp_flag_day_omeq0 == 0)]

    #operation mode: Ocean Glint
    tmp_xco2_day_omeq1 = np.array(in_oco2_nc_data.variables['xco2'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:]==1)]
    tmp_lon_day_omeq1 = np.array(in_oco2_nc_data.variables['longitude'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]
    tmp_lat_day_omeq1 = np.array(in_oco2_nc_data.variables['latitude'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]
    tmp_co2_prior_day_omeq1 = np.array(in_oco2_nc_data.variables['co2_profile_apriori'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1), :][0]
    tmp_xco2_prior_day_omeq1 = np.array(in_oco2_nc_data.variables['xco2_apriori'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]
    tmp_pres_day_omeq1 = np.array(in_oco2_nc_data.variables['pressure_levels'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1), :][0]
    tmp_surfp_day_omeq1 = np.array(in_oco2_nc_data.groups['Retrieval']['psurf'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]
    tmp_xco2_uncert_day_omeq1 = np.array(in_oco2_nc_data.variables['xco2_uncertainty'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]
    tmp_pressure_weight_day_omeq1 = np.array(in_oco2_nc_data.variables['pressure_weight'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1), :][0]
    tmp_xco2_avk_uncert_day_omeq1 = np.array(in_oco2_nc_data.variables['xco2_averaging_kernel'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1), :][0]
    tmp_time_day_omeq1 = np.array(in_oco2_nc_data.variables['time'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]
    tmp_flag_day_omeq1 = np.array(in_oco2_nc_data.variables['xco2_quality_flag'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]
    tmp_lo_flag_omeq1 = np.array(in_oco2_nc_data.groups['Sounding']['land_water_indicator'][:])[np.where(in_oco2_nc_data.groups['Sounding']['operation_mode'][:] == 1)]


    tmp_xco2_day_omeq1_lo_gd = tmp_xco2_day_omeq1[np.where(tmp_flag_day_omeq1==0)]
    tmp_lon_day_omeq1_lo_gd = tmp_lon_day_omeq1[np.where(tmp_flag_day_omeq1 == 0)]
    tmp_lat_day_omeq1_lo_gd = tmp_lat_day_omeq1[np.where(tmp_flag_day_omeq1 == 0)]
    tmp_co2_prior_day_omeq1_lo_gd = tmp_co2_prior_day_omeq1[np.where(tmp_flag_day_omeq1 == 0), :][0]
    tmp_xco2_prior_day_omeq1_lo_gd = tmp_xco2_prior_day_omeq1[np.where(tmp_flag_day_omeq1 == 0)]
    tmp_pres_day_omeq1_lo_gd = tmp_pres_day_omeq1[np.where(tmp_flag_day_omeq1 == 0), :][0]
    tmp_surfp_day_omeq1_lo_gd = tmp_surfp_day_omeq1[np.where(tmp_flag_day_omeq1 == 0)]
    tmp_xco2_uncert_day_omeq1_lo_gd = tmp_xco2_uncert_day_omeq1[np.where(tmp_flag_day_omeq1 == 0)]
    tmp_pressure_weight_day_omeq1_lo_gd = tmp_pressure_weight_day_omeq1[np.where(tmp_flag_day_omeq1 == 0), :][0]
    tmp_xco2_avk_uncert_day_omeq1_lo_gd = tmp_xco2_avk_uncert_day_omeq1[np.where(tmp_flag_day_omeq1 == 0), :][0]
    tmp_time_day_omeq1_lo_gd = tmp_time_day_omeq1[np.where(tmp_flag_day_omeq1 == 0)]
    tmp_lo_flag_omeq1_lo_gd = tmp_lo_flag_omeq1[np.where(tmp_flag_day_omeq1 == 0)]

    #land
    tmp_xco2_day_omeq1_gd = tmp_xco2_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd==0)]
    tmp_lon_day_omeq1_gd = tmp_lon_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0)]
    tmp_lat_day_omeq1_gd = tmp_lat_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0)]
    tmp_co2_prior_day_omeq1_gd = tmp_co2_prior_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0), :][0]
    tmp_xco2_prior_day_omeq1_gd = tmp_xco2_prior_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0)]
    tmp_pres_day_omeq1_gd = tmp_pres_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0), :][0]
    tmp_surfp_day_omeq1_gd = tmp_surfp_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0)]
    tmp_xco2_uncert_day_omeq1_gd = tmp_xco2_uncert_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0)]
    tmp_pressure_weight_day_omeq1_gd = tmp_pressure_weight_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0), :][0]
    tmp_xco2_avk_uncert_day_omeq1_gd = tmp_xco2_avk_uncert_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0), :][0]
    tmp_time_day_omeq1_gd = tmp_time_day_omeq1_lo_gd[np.where(tmp_lo_flag_omeq1_lo_gd == 0)]



    xco2_day = np.concatenate((tmp_xco2_day_omeq0_gd, tmp_xco2_day_omeq1_gd), axis=0)
    lon_day = np.concatenate((tmp_lon_day_omeq0_gd, tmp_lon_day_omeq1_gd), axis=0)
    lat_day = np.concatenate((tmp_lat_day_omeq0_gd, tmp_lat_day_omeq1_gd), axis=0)
    co2_prior_day = np.concatenate((tmp_co2_prior_day_omeq0_gd, tmp_co2_prior_day_omeq1_gd), axis=0)
    xco2_prior_day = np.concatenate((tmp_xco2_prior_day_omeq0_gd, tmp_xco2_prior_day_omeq1_gd), axis=0)
    pres_day = np.concatenate((tmp_pres_day_omeq0_gd, tmp_pres_day_omeq1_gd), axis=0)
    surfp_day = np.concatenate((tmp_surfp_day_omeq0_gd, tmp_surfp_day_omeq1_gd), axis=0)
    xco2_uncert_day = np.concatenate((tmp_xco2_uncert_day_omeq0_gd, tmp_xco2_uncert_day_omeq1_gd), axis=0)
    pressure_weight_day = np.concatenate((tmp_pressure_weight_day_omeq0_gd, tmp_pressure_weight_day_omeq1_gd), axis=0)
    xco2_avk_day = np.concatenate((tmp_xco2_avk_uncert_day_omeq0_gd, tmp_xco2_avk_uncert_day_omeq1_gd), axis=0)
    time_day = np.concatenate((tmp_time_day_omeq0_gd, tmp_time_day_omeq1_gd), axis=0)

    #sort from small to large
    xco2_day = xco2_day[np.argsort(time_day)]
    lon_day = lon_day[np.argsort(time_day)]
    lat_day = lat_day[np.argsort(time_day)]
    co2_prior_day = co2_prior_day[np.argsort(time_day), :]
    xco2_prior_day = xco2_prior_day[np.argsort(time_day)]
    pres_day = pres_day[np.argsort(time_day)]
    surfp_day = surfp_day[np.argsort(time_day)]
    xco2_uncert_day = xco2_uncert_day[np.argsort(time_day)]
    pressure_weight_day = pressure_weight_day[np.argsort(time_day)]
    xco2_avk_day = xco2_avk_day[np.argsort(time_day)]
    time_day = time_day[np.argsort(time_day)]

    #1s average
    index_length = time_day.shape[0]

    # 1s average from good quality data
    tmp_1s_xco2 = np.zeros(index_length)
    tmp_1s_xco2_var = np.zeros(index_length)
    tmp_1s_lat = np.zeros(index_length)
    tmp_1s_lon = np.zeros(index_length)
    tmp_1s_co2_apri_profile = np.zeros((index_length, 20))
    tmp_1s_xco2_apri = np.zeros(index_length)
    tmp_1s_xco2_avk = np.zeros((index_length, 20))
    tmp_1s_pres = np.zeros((index_length, 20))
    tmp_1s_surf = np.zeros(index_length)
    tmp_1s_xco2_pwf = np.zeros((index_length, 20))
    tmp_1s_time = np.zeros(index_length)

    num = 0
    start_index = 0
    if index_length == 0:
        oco2_cur_time_date += datetime.timedelta(days=1)
        continue

    while (True):

        begin_time_num = time_day[start_index]
        # find time interval index in 1s
        I_flag = True  # i exit normally
        for i in range(start_index, index_length):

            time_num = time_day[i]

            if (time_num - begin_time_num) <= 1:
                continue
            else:
                if i == index_length - 1:
                    I_flag = False
                break

        last_index = i
        if i == index_length - 1 and I_flag:
            last_index = i + 1

        tmp_1s_xco2[num], tmp_1s_xco2_var[num], tmp_1s_lat[num], tmp_1s_lon[num], tmp_1s_co2_apri_profile[num, :], tmp_1s_xco2_apri[num], \
        tmp_1s_xco2_avk[num, :], tmp_1s_pres[num, :], tmp_1s_surf[num], tmp_1s_xco2_pwf[num, :], tmp_1s_time[num] \
            = oco2_weighting_calculation(xco2_day, xco2_uncert_day, lat_day, lon_day, co2_prior_day, xco2_prior_day,xco2_avk_day, pres_day, surfp_day, \
                      pressure_weight_day, time_day, start_index, last_index, is_10s=False)

        num += 1

        if i == index_length - 1 and I_flag:
            break

        start_index = i


    # lon_ht_w50 = np.array(tmp_1s_lon)[np.where(np.array(tmp_1s_lon)[:] > -50)]
    # lon_lt_w30 = lon_ht_w50[np.where(lon_ht_w50 < -30)]
    # lat_ht_w50 = np.array(tmp_1s_lat)[np.where(np.array(tmp_1s_lon) > -50)]
    # lat_lt_w30 = lat_ht_w50[np.where(lon_ht_w50 < -30)]
    # lat_ht_n10 = lon_ht_w50[np.where(lat_lt_w30 > 10)]
    #
    # print(lat_ht_n10)
    # exit()
    # 10s calculation
    tmp_10s_xco2, tmp_10s_xco2_var, tmp_10s_lat, tmp_10s_lon, tmp_10s_co2_apri_profile, tmp_10s_xco2_apri, \
    tmp_10s_xco2_avk, tmp_10s_pres, tmp_10s_surf, tmp_10s_xco2_pwf, tmp_10s_time, num_10s \
        = oco2_10s_all(tmp_1s_xco2, tmp_1s_xco2_var, tmp_1s_lat, tmp_1s_lon, \
                       tmp_1s_co2_apri_profile, tmp_1s_xco2_apri, tmp_1s_xco2_avk, tmp_1s_pres, tmp_1s_surf, \
                       tmp_1s_xco2_pwf, tmp_1s_time, num)
    # output


    out_file_path = out_oco2_dir+oco2_filename
    print(out_file_path+' is writing')

    out_oco2 = Dataset(out_file_path, 'w')
    Dataset.createDimension(out_oco2, dimname='Time', size=None)
    Dataset.createDimension(out_oco2, dimname='level', size=20)

    Dataset.createDimension(out_oco2, dimname='DateStrLen', size=23)
    Times_out = Dataset.createVariable(out_oco2, 'Times', datatype=np.float64, dimensions=('Time'))
    Times_out.units = "seconds since 1970-01-01 00:00:00"
    Times_out[:] = tmp_10s_time[:num_10s]


    level_array = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    level_out = Dataset.createVariable(out_oco2, 'level', datatype=np.float32, dimensions=('level'))
    level_out[:] = np.arange(1, 21, 1)

    xco2_out = Dataset.createVariable(out_oco2, 'xco2', datatype=np.float32, dimensions=('Time'))
    xco2_out[:] = tmp_10s_xco2[:num_10s]*10**-6

    # xco2 prior
    xco2_apri_out = Dataset.createVariable(out_oco2, 'xco2_apriori', datatype=np.float32, dimensions=('Time'))
    xco2_apri_out[:] = tmp_10s_xco2_apri[:num_10s]*10**-6

    xco2_uncert_out = Dataset.createVariable(out_oco2, 'xco2_uncert', datatype=np.float32, dimensions=('Time'))
    xco2_uncert_out[:] = tmp_10s_xco2_var[:num_10s]*10**-6

    # lat
    lat_out = Dataset.createVariable(out_oco2, 'latitude', datatype=np.float32, dimensions=('Time'))
    lat_out[:] = tmp_10s_lat[:num_10s]

    lon_out = Dataset.createVariable(out_oco2, 'longitude', datatype=np.float32, dimensions=('Time'))
    lon_out[:] = tmp_10s_lon[:num_10s]

    # pressure
    pressure_out = Dataset.createVariable(out_oco2, 'pressure', datatype=np.float32, dimensions=('Time', 'level'))
    pressure_out[:] = tmp_10s_pres[:num_10s, :]*100


    sup_pressure_out = Dataset.createVariable(out_oco2, 'surface_pressure', datatype=np.float32,dimensions=('Time'))
    sup_pressure_out[:] = tmp_10s_surf[:num_10s]*100

    co2_apri_out = Dataset.createVariable(out_oco2, 'co2_profile_apriori', datatype=np.float32,dimensions=('Time', 'level'))
    co2_apri_out[:] = tmp_10s_co2_apri_profile[:num_10s, :]*10**-6

    xco2_avk_out = Dataset.createVariable(out_oco2, 'xco2_averaging_kernel_matrix', datatype=np.float32,dimensions=('Time', 'level'))
    xco2_avk_out[:] = tmp_10s_xco2_avk[:num_10s, :]

    xco2_pwf_out = Dataset.createVariable(out_oco2, 'xco2_pressure_weighting_function', datatype=np.float32,dimensions=('Time', 'level'))
    xco2_pwf_out[:] = tmp_10s_xco2_pwf[:num_10s, :]

    oco2_cur_time_date += datetime.timedelta(days=1)



pro_end_time = time.time()

print('total running time is:', pro_end_time - pro_begin_time)



