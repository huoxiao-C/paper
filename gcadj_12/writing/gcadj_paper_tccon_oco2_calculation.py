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


tccon = tccon_oco2_pd.iloc[:, 13*5+2][~np.isnan(tccon_oco2_pd.iloc[:, 13*5+2])]*10**6
prior = tccon_oco2_pd.iloc[:, 13*5+3][~np.isnan(tccon_oco2_pd.iloc[:, 13*5+3])]*10**6
posterior = tccon_oco2_pd.iloc[:, 13*5+4][~np.isnan(tccon_oco2_pd.iloc[:, 13*5+4])]*10**6
oco = tccon_oco2_pd.iloc[:, 13*5+1][~np.isnan(tccon_oco2_pd.iloc[:, 13*5+1])]*10**6

print('13: ', np.mean(np.abs(oco-tccon)), np.mean(np.abs(prior-tccon)), np.mean(np.abs(posterior-tccon)))
print('13: ', np.mean(oco-tccon), np.mean(prior-tccon), np.mean(posterior-tccon))
# print('13 oco: ', oco.max(), oco.min(), oco.median(), oco.shape[0], oco)
# print('\noco:', oco, '\n')
# print('\nposterior:', posterior, '\n')
# print('\nprior:', prior, '\n')

tccon = tccon_oco2_pd.iloc[:, 21*5+2][~np.isnan(tccon_oco2_pd.iloc[:, 21*5+2])]*10**6
prior = tccon_oco2_pd.iloc[:, 21*5+3][~np.isnan(tccon_oco2_pd.iloc[:, 21*5+3])]*10**6
posterior = tccon_oco2_pd.iloc[:, 21*5+4][~np.isnan(tccon_oco2_pd.iloc[:, 21*5+4])]*10**6
oco = tccon_oco2_pd.iloc[:, 21*5+1][~np.isnan(tccon_oco2_pd.iloc[:, 21*5+1])]*10**6

# print('21: ', np.mean(tccon), np.mean(prior), np.mean(posterior), np.mean(oco))
# print('21 oco: ', oco.max(), oco.min(), oco.median(), oco.shape[0])
print('21: ', np.mean(np.abs(oco-tccon)), np.mean(np.abs(prior-tccon)), np.mean(np.abs(posterior-tccon)))
print('21: ', np.mean(oco-tccon), np.mean(prior-tccon), np.mean(posterior-tccon))
tccon = tccon_oco2_pd.iloc[:, 5*5+2][~np.isnan(tccon_oco2_pd.iloc[:, 5*5+2])]*10**6
prior = tccon_oco2_pd.iloc[:, 5*5+3][~np.isnan(tccon_oco2_pd.iloc[:, 5*5+3])]*10**6
posterior = tccon_oco2_pd.iloc[:, 5*5+4][~np.isnan(tccon_oco2_pd.iloc[:, 5*5+4])]*10**6
oco = tccon_oco2_pd.iloc[:, 5*5+1][~np.isnan(tccon_oco2_pd.iloc[:, 5*5+1])]*10**6

# print('5: ', np.mean(tccon), np.mean(prior), np.mean(posterior), np.mean(oco))
# print('5 oco: ', oco.max(), oco.min(), oco.median(), oco.shape[0])
print('5: ', np.mean(np.abs(oco-tccon)), np.mean(np.abs(prior-tccon)), np.mean(np.abs(posterior-tccon)))
print('5: ', np.mean(oco-tccon), np.mean(prior-tccon), np.mean(posterior-tccon))