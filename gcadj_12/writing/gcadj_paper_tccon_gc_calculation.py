import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

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



tccon_pd_num = int(tccon_pd.shape[1]/4)
prior_xco2_list = []
posterior_xco2_list = []
tccon_xco2_list = []
for col in range(tccon_pd_num):
    prior_xco2_list.append(tccon_pd.iloc[:, col*4+1][~np.isnan(tccon_pd.iloc[:, col*4+1])])
    posterior_xco2_list.append(tccon_pd.iloc[:, col*4+2][~np.isnan(tccon_pd.iloc[:, col*4+2])])
    tccon_xco2_list.append(tccon_pd.iloc[:, col*4+3][~np.isnan(tccon_pd.iloc[:, col*4+3])])

# plt.scatter(range(len(prior_xco2_list[6])), prior_xco2_list[6])
# plt.scatter(range(len(prior_xco2_list[6])), posterior_xco2_list[6])
# plt.scatter(range(len(prior_xco2_list[6])), tccon_xco2_list[6])
# plt.show()