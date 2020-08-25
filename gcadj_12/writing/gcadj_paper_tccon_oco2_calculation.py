import pandas as pd
import numpy as np



tccon_oco2_pd = pd.read_excel('/home/huoxiao/n1s1/huoxiao_share/paper_data/tccon_oco2_nd.xlsx', index_col=0)
tccon_oco2_all_name_list = ['Sites/Bialystok', 'Bremen', 'Burgos', 'Caltech',
                            'Darwin', 'Dryden', 'East Trout Lake', 'Garmisch',
                            'Izana', 'Pasadena', 'Saga', 'Karlsruhe',
                            'Lamont', 'Orleans', 'Park Falls', 'Rikubetsu',
                            'Tsukuba', 'Wollongong', 'Zugspitze']
for i in range(30):
    print(tccon_oco2_pd.columns[i*3])
tccon_oco2_pd_num = int(tccon_oco2_pd.shape[1]/3)
oco2_xco2 = tccon_oco2_pd.iloc[:, 3*18+1][~np.isnan(tccon_oco2_pd.iloc[:, 3*18+1])]
tccon_xco2 = tccon_oco2_pd.iloc[:, 3*18+2][~np.isnan(tccon_oco2_pd.iloc[:, 3*18+2])]