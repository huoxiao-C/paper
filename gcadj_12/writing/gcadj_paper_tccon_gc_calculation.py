import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.geoaxes import GeoAxes
from matplotlib.cm import get_cmap
from matplotlib.ticker import MaxNLocator
import matplotlib.dates as mdates
import datetime
import matplotlib.ticker as ticker

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

class TCCON_GC():

    def __init__(self, tccon_site_name, time,
                 pri_xco2, pos_xco2, tccon_xco2):
        self.time = time
        self.tccon_site_name = tccon_site_name
        self.pri_xco2 = pri_xco2
        self.pos_xco2 = pos_xco2
        self.tccon_xco2 = tccon_xco2

    def tccon_gc_plot(self, tccon_site_plot_index,
                      tick_font_size,
                      label_font_size):
        fig, axes = plt.subplots(5, 3, figsize=(32, 30))
        for row in range(5):
            for col in range(3):
                print(row, col, row*3+col)
                ax = axes[row, col]
                ax.scatter(self.time, self.pri_xco2[tccon_site_plot_index[row*3+col]], label='Prior')
                ax.scatter(self.time, self.pos_xco2[tccon_site_plot_index[row*3+col]], label='Posterior')
                ax.scatter(self.time, self.tccon_xco2[tccon_site_plot_index[row*3+col]], label='TCCON')

                formatter = mdates.DateFormatter('%Y-%m')



                # ax.xaxis.set_major_formatter(ymsFmt)
                # ax.xaxis.set_minor_formatter(months)
                ax.xaxis.set_major_formatter(formatter)
                daylocator = mdates.DayLocator()
                monthlocator = mdates.MonthLocator([i for i in range(1, 13, 2)])
                ax.xaxis.set_major_locator(monthlocator)
                ax.xaxis.set_minor_locator(daylocator)
                ax.set_title(self.tccon_site_name[tccon_plot_index_list[row*3+col]], fontsize=label_font_size)
                # ax.set_xticks()
                ax.set_yticks(np.arange(397, 418, 2))

                if row==0 and col ==0:
                    ax.legend(fontsize=tick_font_size, loc='lower left', prop={'size': label_font_size})
        plt.savefig('/home/huoxiao/paper_data/tccon.svg', dpi=300, format='svg', bbox_inches = 'tight', pad_inches = 0.5)


    def tccon_gc_box_plot(self, tccon_site_plot_index,
                      tick_font_size,
                      label_font_size):

        prior_xco2_list = []
        posterior_xco2_list = []
        tccon_xco2_list = []
        tccon_site_name_list = []
        colors = ['pink', 'lightblue', 'lightgreen']
        for site_index in tccon_site_plot_index:
            prior_xco2_list.append(self.pri_xco2[site_index][~np.isnan(self.pri_xco2[site_index])])
            posterior_xco2_list.append(np.array(self.pos_xco2[site_index][~np.isnan(self.pos_xco2[site_index])]))
            tccon_xco2_list.append(self.tccon_xco2[site_index][~np.isnan(self.tccon_xco2[site_index])])
            tccon_site_name_list.append(self.tccon_site_name[site_index])

        fig, axes = plt.subplots(1, 1, figsize=(32, 16))
        pri_bplot = axes.boxplot(prior_xco2_list, positions=np.arange(2, 92, 6),
                                 widths=1, patch_artist=True, boxprops=dict(facecolor=colors[0]))
        pos_bplot = axes.boxplot(posterior_xco2_list, positions=np.arange(3, 93, 6),
                                 widths=1, patch_artist=True, boxprops=dict(facecolor=colors[1]))
        tccon_bplot = axes.boxplot(tccon_xco2_list, positions=np.arange(4, 94, 6),
                                   widths=1, patch_artist=True, boxprops=dict(facecolor=colors[2]))

        axes.set_xticks(np.arange(3, 93, 6))
        axes.set_yticks(np.arange(395, 418, 2))
        axes.tick_params(axis="y", labelsize=tick_font_size)
        axes.set_xticklabels(tccon_site_name_list, fontsize=tick_font_size)
        axes.set_xlim(1, 89)
        axes.margins(x=100)
        axes.legend([pri_bplot["boxes"][0], pos_bplot["boxes"][0], tccon_bplot["boxes"][0]],
                    ['Prior', 'Posterior', 'TCCON'], loc='upper right', fontsize=tick_font_size)
        plt.savefig('/home/huoxiao/paper_data/tccon_boxplot.svg', dpi=300, format='svg', bbox_inches='tight', pad_inches=0.5)
        # for pri_patch,pos_patch, tccon_patch in \
        #         zip(pri_bplot['boxes'],
        #         pos_bplot['boxes'],
        #         tccon_bplot['boxes']):
        #     pri_patch.set_facecolor(colors[0])
        #     pos_patch.set_facecolor(colors[1])
        #     tccon_patch.set_facecolor(colors[2])
tccon_pd = pd.read_excel('/home/huoxiao/n1s1/huoxiao_share/paper_data/tccon_gc_single_site_nd.xlsx', index_col=0)
#oco 405.5815473130236 #tccon 406.32164634119437
#the Southern Hemisphere: Ascension,Darwin, Lauder, Reunion, Wollongong
#north pole:Eureka,  Ny Alesund
tccon_all_name_list = ['Ascension', 'Anmyeondo', 'Sites/Bialystok', 'Bremen',
                       'Burgos', 'Caltech', 'Darwin', 'Dryden', 'East Trout Lake',
                       'Eureka', 'Garmisch', 'Izana', 'Pasadena', 'Saga', 'Karlsruhe',
                       'Lauder',  'Lauder', 'Lamont', 'Orleans', 'Park Falls', 'Paris',
                       'Reunion', 'Rikubetsu', 'Sodankyla', 'Ny Alesund', 'Tsukuba',
                       'Wollongong', 'Zugspitze']


#13 5 21
begin_time_str = '2018-01-01 00:00:00'
begin_time = datetime.datetime.strptime(begin_time_str, '%Y-%m-%d %H:%M:%S')
end_time = datetime.datetime.strptime('2019-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')
tccon_pd_num = int(tccon_pd.shape[1]/4)

prior_xco2_list = []
posterior_xco2_list = []
tccon_xco2_list = []
for col in range(tccon_pd_num):
    prior_site_array = np.full_like(np.arange(365*24), np.nan, dtype=np.double)
    posterior_site_array = np.full_like(np.arange(365 * 24), np.nan, dtype=np.double)
    tccon_site_array = np.full_like(np.arange(365 * 24), np.nan, dtype=np.double)
    test_num = 0
    for site_num in range(tccon_pd.iloc[:, col*4+1].shape[0]):
        if np.isnan(tccon_pd.iloc[site_num, col*4+1]):
            continue
        site_date = datetime.datetime.strptime(tccon_pd.iloc[site_num, col*4+0], '%Y-%m-%d %H:%M:%S')
        hour_index = int((site_date-begin_time).total_seconds()/3600)
        test_num +=1

        prior_site_array[hour_index] = tccon_pd.iloc[site_num, col*4+1]*10**6
        posterior_site_array[hour_index] = tccon_pd.iloc[site_num, col * 4 + 2]*10**6
        tccon_site_array[hour_index] = tccon_pd.iloc[site_num, col * 4 + 3] * 10 ** 6
    prior_xco2_list.append(prior_site_array)
    posterior_xco2_list.append(posterior_site_array)
    tccon_xco2_list.append(tccon_site_array)

tccon_plot_index_list = [3, 6, 8, 10, 13, 14, 15, 17, 18, 19, 22, 23, 25, 26, 27]
time_list = [begin_time +datetime.timedelta(hours=i) for i in range(365*24)]
#plot
tick_font_size = 18
label_font_size = 30

TCCON2GC = TCCON_GC(tccon_all_name_list, time_list,
                    prior_xco2_list, posterior_xco2_list,
                    tccon_xco2_list)
# TCCON2GC.tccon_gc_plot(tccon_plot_index_list, tick_font_size, label_font_size)
TCCON2GC.tccon_gc_box_plot(tccon_plot_index_list, tick_font_size, label_font_size)

plt.show()

# plt.scatter(range(len(prior_xco2_list[6])), prior_xco2_list[6])
# plt.scatter(range(len(prior_xco2_list[6])), posterior_xco2_list[6])
# plt.scatter(range(len(prior_xco2_list[6])), tccon_xco2_list[6])
# plt.show()
