from netCDF4 import Dataset, num2date, date2num
import numpy as np
import datetime
import os
import subprocess
import ctypes
import time
import calendar
import datetime


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


def pressure_weigt_function(pres, surfp, level, lat_num, lon_num):
    ''''

    '''
    h = np.zeros((level, lat_num, lon_num))
    for i in range(level):
        if i == 0:
            h[i, :, :] = np.abs(-pres[i, :, :] + (pres[i + 1, :, :] - pres[i, :, :]) / np.log(
                pres[i + 1, :, :] / pres[i, :, :])) / surfp[:, :]
        elif i == level - 1:
            h[i, :, :] = np.abs(pres[i, :, :] - (pres[i, :, :] - pres[i - 1, :, :]) / np.log(
                pres[i, :, :] / pres[i - 1, :, :])) / surfp[:, :]
        else:
            h[i, :, :] = np.abs(-pres[i, :, :] + (pres[i + 1, :, :] - pres[i, :, :]) / np.log(pres[i + 1, :, :] \
                                                                                              / pres[i, :, :]) + pres[i,
                                                                                                                 :,
                                                                                                                 :] - (
                                            pres[i, :, :] - pres[i - 1, :, :]) / np.log(
                pres[i, :, :] / pres[i - 1, :, :])) / surfp[:, :]

    return h


def geos_chem_grid_axis():
    ''''
    calculate geos_chem grid longitude center and latitude center
    '''

    lon = np.zeros(144)
    lat = np.zeros(92)
    for i in range(144):
        lon[i] = 2.5 * i - 180.0

    for j in range(91):
        lat[j] = 2 * j - 90.0

    return lon, lat


def add_element(data, time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, xco2_pri_list, \
                xco2_avk_list, pre_list, sufp_list, xco2_pwf_list, j):
    '''
    add element to list
    :param data:
    :param time_list:
    :param co2_list:
    :param pre_list:
    :param lat_list:
    :param lon_list:
    :param pri_list:
    :param avk_list:
    :param Sinv_list:
    :return:
    '''
    # try:
    #     time_list.append(data.variable['Times'][j])
    # except:
    #     print(j)
    #     print(data.variables['Times'][0])
    #     exit(1)
    if data.variables['longitude'][j] > -50 and data.variables['longitude'][j] < -30 and \
            data.variables['latitude'][j] > 10:
        print('lon, lat', data.variables['latitude'][j], data.variables['latitude'][j])

    time_list.append(data.variables['Times'][j])

    pre_list.append(data.variables['pressure'][j, :])

    co2_pri_list.append(data.variables['co2_profile_apriori'][j, :])

    xco2_list.append(data.variables['xco2'][j])

    lat_list.append(data.variables['latitude'][j])

    lon_list.append(data.variables['longitude'][j])

    xco2_avk_list.append(data.variables['xco2_averaging_kernel_matrix'][j, :])

    sufp_list.append((data.variables['surface_pressure'][j]))

    xco2_uncert_list.append(data.variables['xco2_uncert'][j])

    xco2_pri_list.append(data.variables['xco2_apriori'][j])

    xco2_avk_list.append(data.variables['xco2_averaging_kernel_matrix'][j, :])
    xco2_pwf_list.append(data.variables['xco2_pressure_weighting_function'][j, :])

    return time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, xco2_pri_list, xco2_avk_list, \
           pre_list, sufp_list, xco2_pwf_list


def extractlist(xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, \
                xco2_pri_list, xco2_avk_list, pre_list, sufp_list, xco2_pwf_list):
    length = len(lon_list)
    pressures = np.zeros((length, 20))
    lats = np.zeros(length)
    lons = np.zeros(length)
    co2_priors = np.zeros((length, 20))
    sufps = np.zeros(length)
    xco2s = np.zeros(length)
    xco2_uncerts = np.zeros(length)
    xco2_priors = np.zeros(length)
    xco2_avks = np.zeros((length, 20))

    xco2_pwf = np.zeros((length, 20))
    for i in range(length):
        xco2s[i] = xco2_list[i]
        pressures[i, :] = pre_list[i]
        lats[i] = lat_list[i]
        lons[i] = lon_list[i]
        co2_priors[i, :] = co2_pri_list[i]
        sufps[i] = sufp_list[i]
        xco2_uncerts[i] = xco2_uncert_list[i]
        xco2_priors[i] = xco2_pri_list[i]
        xco2_avks[i, :] = xco2_avk_list[i]
        xco2_pwf[i, :] = xco2_pwf_list[i]

    return xco2s, xco2_uncerts, lats, lons, co2_priors, xco2_priors, xco2_avks, pressures, sufps, xco2_pwf


def get_ij(lon, lat, delta_x, delta_y, NX, NY):
    ''''
    '''
    tlon = int((lon + 180.0) / delta_x + 1.5)

    tlat = int((lat + 90.0) / delta_y + 1.5)

    if (tlon > NX):
        tlon = tlon - NX

    return tlon, tlat


def tansat_grid_to_geoschem(time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, \
                            xco2_pri_list, xco2_avk_list, pre_list, sufp_list, xco2_pwf_list, delta_x, delta_y, NX, NY):
    ''''
    '''

    tmp_time_list = []
    tmp_pre_list = []
    tmp_sufp_list = []
    tmp_xco2_list = []
    tmp_lat_list = []
    tmp_lon_list = []
    tmp_co2_pri_list = []
    tmp_xco2_pwf_list = []
    tmp_xco2_pri_list = []
    tmp_xco2_avk_list = []
    tmp_xco2_uncert_list = []

    ij_list = []

    gc_lon, gc_lat = geos_chem_grid_axis()

    for index in range(len(time_list)):
        lon = lon_list[index]
        lat = lat_list[index]
        i, j = get_ij(lon, lat, delta_x, delta_y, NX, NY)
        ij_list.append('{:0>3d}+{:0>3d}'.format(i, j))

    ij_set = set(ij_list)
    # average on the sanme grid
    for index in range(len(ij_set)):
        subs_list = [i for i, item in enumerate(ij_list) if item == list(ij_set)[index]]

        #
        list_len = len(subs_list)
        tmp_time = 0.0
        tmp_pre = np.zeros(20)
        tmp_sufp = 0.0
        tmp_xco2 = 0.0
        tmp_lat = 0.0
        tmp_lon = 0.0
        tmp_xco2_pri = 0.0
        tmp_xco2_uncert = 0.0
        tmp_xco2_pwf = np.zeros(20)
        tmp_xco2_avk = np.zeros(20)
        tmp_co2_pri = np.zeros(20)

        for index1 in range(list_len):
            list_index = subs_list[index1]
            tmp_time += time_list[list_index] / list_len
            tmp_pre += pre_list[list_index] / list_len
            tmp_sufp += sufp_list[list_index] / list_len
            tmp_xco2 += xco2_list[list_index] / list_len
            tmp_xco2_pri += xco2_pri_list[list_index] / list_len
            tmp_xco2_uncert += xco2_uncert_list[list_index] / list_len
            tmp_xco2_avk += xco2_avk_list[list_index] / list_len
            tmp_xco2_pwf += xco2_pwf_list[list_index] / list_len
            # tmp_lat += lat_list[list_index]/list_len
            # tmp_lon += lon_list[list_index]/list_len
            tmp_co2_pri += co2_pri_list[list_index] / list_len

        i = int(list(ij_set)[index][0:3])

        j = int(list(ij_set)[index][4:7])

        tmp_lon = gc_lon[i - 1]
        tmp_lat = gc_lat[j - 1]

        tmp_time_list.append(tmp_time)
        tmp_pre_list.append(tmp_pre)
        tmp_sufp_list.append(tmp_sufp)
        tmp_xco2_list.append(tmp_xco2)
        tmp_lat_list.append(tmp_lat)
        tmp_lon_list.append(tmp_lon)
        tmp_co2_pri_list.append(tmp_co2_pri)
        tmp_xco2_pri_list.append(tmp_xco2_pri)
        tmp_xco2_uncert_list.append(tmp_xco2_uncert)
        tmp_xco2_avk_list.append(tmp_xco2_avk)
        tmp_xco2_pwf_list.append(tmp_xco2_pwf)
    return tmp_time_list, tmp_xco2_list, tmp_xco2_uncert_list, tmp_lat_list, tmp_lon_list, tmp_co2_pri_list, \
           tmp_xco2_pri_list, tmp_xco2_avk_list, tmp_pre_list, tmp_sufp_list, tmp_xco2_pwf_list


def observation_error_variance_bias_mean(begin_time, end_time, oco2_dir, gc_dir, gc_level, delta_x, delta_y, NX, NY):
    '''

    :param begin_time:
    :param end_time:
    :param inoco2_dir:
    :param gc_level:
    :param delta_x:
    :param delta_y:
    :param NX:
    :param NY:
    :return:
    '''

    observation_num = np.zeros((NY, NX))
    ym_minus_yo = np.zeros((NY, NX))
    ym_minus_yo_bar = np.zeros((NY, NX))
    yo = np.zeros((NY, NX))
    yo_bar = np.zeros((NY, NX))
    # mean difference between model observation and satellite observation

    while (end_time != begin_time):
        print('observation error variance stage1', datetime.datetime.strftime(begin_time, '%Y%m%d%H%M%S'))

        oco2_file_path = oco2_dir + 'oco2_' + begin_time.strftime(
            '%Y%m%d%H%M%S') + '.nc'

        if not subprocess.call('test -f ' + oco2_file_path, shell=True):

            oco2_nc_data = Dataset(oco2_file_path, 'r')

            tmp_xco2_file = oco2_nc_data.variables['xco2'][:]
            tmp_lat_file = oco2_nc_data.variables['latitude'][:]
            tmp_lon_file = oco2_nc_data.variables['longitude'][:]

            gc_spc_filename = 'GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour, begin_time.minute)
            gc_met_filename = 'GEOSChem.StateMet.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
                .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour, begin_time.minute)
            gc_spc_path = gc_dir + gc_spc_filename
            gc_met_path = gc_dir + gc_met_filename
            gc_spc_nc_data = Dataset(gc_spc_path, 'r')
            gc_met_nc_data = Dataset(gc_met_path, 'r')

            # convert profile concentration to column concentration
            h = pressure_weigt_function(gc_met_nc_data.variables['Met_PMID'][0, :, :, :], \
                                        gc_met_nc_data.variables['Met_PMID'][0, 0, :, :], gc_level, NY, NX)
            gc_xco2 = np.sum(h[:, :, :] * (gc_spc_nc_data.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
            # search longitude latitude
            for oco2ith in range(tmp_xco2_file.shape[0]):
                lon_index, lat_index = get_ij(tmp_lon_file[oco2ith], tmp_lat_file[oco2ith], delta_x, delta_y, NX, NY)
                lon_index = lon_index - 1
                lat_index = lat_index - 1
                observation_num[lat_index, lon_index] += 1
                ym_minus_yo[lat_index, lon_index] = ym_minus_yo[lat_index, lon_index] \
                                                    + gc_xco2[lat_index, lon_index] - tmp_xco2_file[oco2ith]
                yo[lat_index, lon_index] += tmp_xco2_file[oco2ith]

        begin_time = begin_time + datetime.timedelta(minutes=20)

    ym_minus_yo_bar[np.where(observation_num != 0)] = ym_minus_yo[np.where(observation_num != 0)] / observation_num[
        np.where(observation_num != 0)]
    yo_bar[np.where(observation_num != 0)] = yo[np.where(observation_num != 0)] / observation_num[
        np.where(observation_num != 0)]
    return ym_minus_yo_bar, yo_bar


def observation_error_variance(begin_time, oco2_dir, gc_dir, gc_level, ym_minus_yo_bar, yo_bar):
    '''

    :param begin_time:
    :param end_time:
    :param inoco2_dir:
    :param gc_level:
    :param ym_minus_yo_bar:
    :param yo_bar:
    :return:
    '''

    xco2_uncert_list = []
    print('observation error variance stage2', datetime.datetime.strftime(begin_time, '%Y%m%d%H%M%S'))
    oco2_file_path = oco2_dir + 'oco2_' + begin_time.strftime('%Y%m%d%H%M%S') + '.nc'
    print('observation error variance stage2', oco2_file_path)
    if not subprocess.call('test -f ' + oco2_file_path, shell=True):

        oco2_nc_data = Dataset(oco2_file_path, 'r')

        tmp_xco2_file = oco2_nc_data.variables['xco2'][:]
        tmp_lat_file = oco2_nc_data.variables['latitude'][:]
        tmp_lon_file = oco2_nc_data.variables['longitude'][:]
        tmp_xco2_uncert_file = oco2_nc_data.variables['xco2_uncert'][:]

        gc_spc_filename = 'GEOSChem.SpeciesConc.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
            .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour, begin_time.minute)
        gc_met_filename = 'GEOSChem.StateMet.{:0>4d}{:0>2d}{:0>2d}_{:0>2d}{:0>2d}z.nc4' \
            .format(begin_time.year, begin_time.month, begin_time.day, begin_time.hour, begin_time.minute)
        gc_spc_path = gc_dir + gc_spc_filename
        gc_met_path = gc_dir + gc_met_filename
        gc_spc_nc_data = Dataset(gc_spc_path, 'r')
        gc_met_nc_data = Dataset(gc_met_path, 'r')

        # convert profile concentration to column concentration
        h = pressure_weigt_function(gc_met_nc_data.variables['Met_PMID'][0, :, :, :], \
                                    gc_met_nc_data.variables['Met_PMID'][0, 0, :, :], gc_level, NY, NX)
        gc_xco2 = np.sum(h[:, :, :] * (gc_spc_nc_data.variables['SpeciesConc_CO2'][0, :, :, :]), axis=0)
        # search longitude latitude
        for oco2ith in range(tmp_xco2_file.shape[0]):
            ym = gc_xco2

            lon_index, lat_index = get_ij(tmp_lon_file[oco2ith], tmp_lat_file[oco2ith], delta_x, delta_y, NX, NY)
            lon_index = lon_index - 1
            lat_index = lat_index - 1
            residual_error = ym[lat_index, lon_index] - tmp_xco2_file[oco2ith] - ym_minus_yo_bar[lat_index, lon_index]
            rrsd = residual_error / yo_bar[lat_index, lon_index]
            observation_error_variance = rrsd * tmp_xco2_file[oco2ith]
            xco2_uncert_list.append(max(abs(observation_error_variance), tmp_xco2_uncert_file[oco2ith]))

    return xco2_uncert_list


def tansat_file_create(filename, time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, \
                       xco2_pri_list, xco2_avk_list, pre_list, sufp_list, xco2_pwf_list):
    ''''
    create oco2 file
    '''

    # exchange column
    exchange_list = [i for i in range(19, -1, -1)]

    output = Dataset(filename, 'w', format='NETCDF3_64BIT_OFFSET')
    Dataset.createDimension(output, dimname='Time', size=None)
    Dataset.createDimension(output, dimname='level', size=20)
    Dataset.createDimension(output, dimname='DateStrLen', size=23)
    Dataset.createDimension(output, dimname='Const', size=1)

    Times_out = Dataset.createVariable(output, 'Times', datatype=np.float64, \
                                       dimensions=('Time'))
    Times_out.units = 'seconds since 1970-01-01 00:00:00'
    Times_out[:] = np.array(time_list)

    xco2s, xco2_uncerts, lats, lons, co2_priors, xco2_priors, xco2_avks, pressures, sufps, xco2_pwf = \
        extractlist(xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, \
                    xco2_pri_list, xco2_avk_list, pre_list, sufp_list, xco2_pwf_list)

    # xco2
    xco2_out = Dataset.createVariable(output, 'xco2', datatype=np.float32, dimensions=('Time'))
    xco2_out[:] = xco2s

    # pressure
    pre_out = Dataset.createVariable(output, 'pressure', datatype=np.float32, dimensions=('Time', 'level'))
    pre_out[:, :] = np.array(pressures)[:, exchange_list]  # column exchange

    # length
    len_out = Dataset.createVariable(output, 'length', datatype=np.float32, dimensions=('Const'))
    len_out[:] = len(lats)
    # latitude
    lat_out = Dataset.createVariable(output, 'latitude', datatype=np.float32, dimensions=('Time'))
    lat_out[:] = lats

    # longitude
    lon_out = Dataset.createVariable(output, 'longitude', datatype=np.float32, dimensions=('Time'))
    lon_out[:] = lons

    # prior profile
    pri_out = Dataset.createVariable(output, 'co2_profile_apriori', datatype=np.float32, dimensions=('Time', 'level'))
    pri_out[:, :] = co2_priors[:, exchange_list]

    # surface pressure
    sufp_out = Dataset.createVariable(output, 'surface_pressure', datatype=np.float32, dimensions=('Time'))
    sufp_out[:] = sufps

    # xco2_uncert
    xco2_uncert_out = Dataset.createVariable(output, 'xco2_uncert', datatype=np.float32,
                                             dimensions=('Time'))
    xco2_uncert_out[:] = xco2_uncerts

    # xco2_prior
    xco2_pri_out = Dataset.createVariable(output, 'xco2_apriori', datatype=np.float32,
                                          dimensions=('Time'))
    xco2_pri_out[:] = xco2_priors

    # xco2 averaging kernel matrix
    xco2_avk_out = Dataset.createVariable(output, 'xco2_avg_kernel', datatype=np.float32,
                                          dimensions=('Time', 'level'))
    xco2_avk_out[:, :] = xco2_avks[:, exchange_list]

    # pressure weighting function
    pwf_out = Dataset.createVariable(output, 'xco2_pressure_weighting_function', datatype=np.float32,
                                     dimensions=('Time', 'level'))
    pwf_out[:, :] = xco2_pwf[:, exchange_list]

    # check
    if np.where(pwf_out[:] > 1)[0].shape[0] != 0 or np.where(pwf_out[:] < -1)[0].shape[0] != 0:
        print('pressure weighting function error')
        print(pwf_out[:, :])
        exit()


begin_year = 2018
begin_month = 1
begin_day = 1
# geos chem grid
delta_x = 2.5
delta_y = 2
NX = 144
NY = 91
gc_level = 47

oco2_time_unit = 'seconds since 1970-01-01 00:00:00'

oco2_begin_time_str = "2017-10-01 00:00:00"
oco2_end_time_str = "2019-01-01 00:00:00"
oco2_begin_time_date = datetime.datetime.strptime(oco2_begin_time_str, '%Y-%m-%d %H:%M:%S')
oco2_end_time_date = datetime.datetime.strptime(oco2_end_time_str, '%Y-%m-%d %H:%M:%S')

gc_file_dir = '/share/nas1_share1/huoxiao_share/GEOSChem_Simulation_test/'
oco210s_dir = '/share/nas1_share1/huoxiao_share/CO2_data/Observation_data/satellite_10s_average/OCO2_Litev9_ND_plus_Glintland/'
oco2pos_dir = '/share/nas1_share1/huoxiao_share/CO2_data/Observation_data/satellite_on_gc_2x25grid/OCO2_CO2_Litev9_ND_plus_GlintLand/'
os.chdir(oco210s_dir)
first_index = 0
oco2_cur_time_date = oco2_begin_time_date

ym_minus_yo_bar, yo_bar = observation_error_variance_bias_mean \
    (oco2_begin_time_date, oco2_end_time_date, oco2pos_dir, gc_file_dir, gc_level, delta_x, delta_y, NX, NY)
# day
while (oco2_cur_time_date != oco2_end_time_date):

    # Determine if the file exists,
    if subprocess.call(
            'test -f ' + oco210s_dir + '{:0>4d}/oco2_LtCO2_{:0>2d}{:0>2d}{:0>2d}'.format(oco2_cur_time_date.year, \
                                                                                         oco2_cur_time_date.year % 2000, \
                                                                                         oco2_cur_time_date.month, \
                                                                                         oco2_cur_time_date.day) + '_B9003r*',
            shell=True):
        oco2_cur_time_date += datetime.timedelta(days=1)
        continue

    command_output = subprocess.check_output(
        "ls  " + '{:0>4d}/oco2_LtCO2_{:0>2d}{:0>2d}{:0>2d}'.format(oco2_cur_time_date.year, \
                                                                   oco2_cur_time_date.year % 2000, \
                                                                   oco2_cur_time_date.month, \
                                                                   oco2_cur_time_date.day) + '_B9003r*', shell=True)

    oco2_filename = files_string2list(command_output.decode())[0]
    oco2_file_path = oco210s_dir + oco2_filename
    print(oco2_file_path)

    oco2_nc_data = Dataset(oco2_file_path, 'r')

    begin_time_file = oco2_cur_time_date

    while (begin_time_file != (oco2_cur_time_date + datetime.timedelta(days=1))):

        t1 = time.time()
        time_list = []
        pre_list = []
        sufp_list = []
        xco2_list = []
        lat_list = []
        lon_list = []
        co2_pri_list = []
        xco2_avk_list = []
        xco2_pwf_list = []
        xco2_pri_list = []
        xco2_uncert_list = []

        time_float_b = oco2_nc_data.variables['Times'][0]
        time_float_e = oco2_nc_data.variables['Times'][-1]

        for num in range(oco2_nc_data.variables['Times'][:].shape[0]):
            time_date = oco2_nc_data.variables['Times'][num]

            # truncate latitude > 60 deg and latitude < -60 deg
            if (oco2_nc_data.variables['latitude'][num] > 60 or oco2_nc_data.variables['latitude'][num] < -60):
                continue

            if ((time_date - date2num(begin_time_file, oco2_time_unit)) / 60 <= 30 and (
                    date2num(begin_time_file, oco2_time_unit) - time_date) / 60 <= 30):
                # print(num2date(time_date, oco2_time_unit), begin_time_file, date2num(begin_time_file, oco2_time_unit),time_date)
                time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, xco2_pri_list, xco2_avk_list, \
                pre_list, sufp_list, xco2_pwf_list \
                    = add_element(oco2_nc_data, time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, \
                                  co2_pri_list, xco2_pri_list, xco2_avk_list, \
                                  pre_list, sufp_list, xco2_pwf_list, num)

        if len(time_list) != 0:
            # oco2 grid to geoschem grid
            time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, co2_pri_list, xco2_pri_list, xco2_avk_list, \
            pre_list, sufp_list, xco2_pwf_list = tansat_grid_to_geoschem(time_list, xco2_list, xco2_uncert_list, \
                                                                         lat_list, lon_list, co2_pri_list, \
                                                                         xco2_pri_list, xco2_avk_list, pre_list, \
                                                                         sufp_list, xco2_pwf_list, \
                                                                         delta_x, delta_y, NX, NY)

            oco2_out_file_path = "/share/nas1_share1/huoxiao_share/CO2_data/Observation_data/satellite_on_gc_2x25grid/" \
                                 + "OCO2_CO2_Litev9_ND_plus_GlintLand_Observation_Error_Covariance/oco2_" + \
                                 datetime.datetime.strftime(begin_time_file, '%Y%m%d%H%M%S') + ".nc"
            print(oco2_out_file_path)
            xco2_uncert_list = observation_error_variance(begin_time_file, oco2pos_dir, gc_file_dir, gc_level, \
                                                          ym_minus_yo_bar, yo_bar)
            tansat_file_create(oco2_out_file_path, time_list, xco2_list, xco2_uncert_list, lat_list, lon_list, \
                               co2_pri_list, \
                               xco2_pri_list, xco2_avk_list, pre_list, sufp_list, xco2_pwf_list)

        begin_time_file = begin_time_file + datetime.timedelta(minutes=60)

    oco2_cur_time_date += datetime.timedelta(days=1)


