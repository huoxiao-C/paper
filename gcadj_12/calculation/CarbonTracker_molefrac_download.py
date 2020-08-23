import subprocess
import datetime




begin_time_str = "2016-10-01 00:00:00"
end_time_str = "2018-01-01 00:00:00"

begin_time_date = datetime.datetime.strptime(begin_time_str, '%Y-%m-%d %H:%M:%S')
end_time_date = datetime.datetime.strptime(end_time_str, '%Y-%m-%d %H:%M:%S')

url_prefix = 'ftp://aftp.cmdl.noaa.gov/products/carbontracker/co2/CT2019/molefractions/co2_total/'
out_dir = '/share/nas1_share1/huoxiao_share/CO2_data/CarbonTracker/CO2_Concentration/2016_2017/ '
while(begin_time_date!=end_time_date):
    url = url_prefix+'CT2019.molefrac_glb3x2_{:0>4d}-{:0>2d}-{:0>2d}.nc'.\
                    format(begin_time_date.year, begin_time_date.month, begin_time_date.day)


    subprocess.call('wget -P ' + out_dir + url, shell=True)
    begin_time_date += datetime.timedelta(days=1)


# #ocean flux
# url_prefix = 'ftp://aftp.cmdl.noaa.gov/products/carbontracker/co2/CT2019/fluxes/three-hourly/'
# out_dir = '~/n1s1/wn12_hx/CO2_data/CarbonTracker/Posterious_flux/ '
# while(begin_time_date!=end_time_date):
#     url = url_prefix+'CT2019.flux1x1.{:0>4d}{:0>2d}{:0>2d}.nc'.\
#                     format(begin_time_date.year, begin_time_date.month, begin_time_date.day)
#
#     print(url)
#     subprocess.call('wget -P ' + out_dir + url, shell=True)
#     begin_time_date += datetime.timedelta(days=1)
#
#
