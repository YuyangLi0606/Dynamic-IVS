#!/usr/bin/env python
# coding: utf-8

# In[216]:


import numpy as np
import math
import pandas as pd
from IPython.core.debugger import set_trace
import multiprocessing
from sklearn.linear_model import HuberRegressor
from sklearn.metrics import mean_squared_error
import sys 
#sys.path.append('C:\\Users\\Bran\\Dropbox\\Research\\Research_Paper\\Implied_vol\\my_codes\\python_code')
import csv 
from datetime import datetime
import pdb
#from skopt.plots import plot_convergence, plot_evaluations, plot_objective

#sys.path.append('/Users/dell/Dropbox/Research/New_data')



 #train_m1, vali_m1, test_m1 = np.split(dat_m1, 3)
#%%  
from nova_test_rolling_scheme_R2_real_alphas import expanding_window
import time 
#from pathos.pools import ProcessPool as Pool 
#import pyspark
#from pathos.pools import ThreadPool as Pool   # Faster for small data, but could be slower for large data ? 
#from pathos.pools import ParallelPool as Pool # Fastest for small data -> can be slow on huge data 


#import tensorflow as tf 
#print(tf.test.gpu_device_name())
# Nested parallel: Not enough space if the data is too big 
start = time.time()
result_pca = []
result_pls = []
result_Enet = []
result_RF = []
result_XGB = []

## the inner dot selection
#K_m = 3
#K_tau = 3
#
#
## fixed degree (AT: degree in te(), degree - 1 setup, i.e. poly_d = 2 represents cubic spline in mgcv te())
#d_m = 3
#d_tau = 2
#
# same for x, z, xz
knots_K = 3 
degree_m = 3


# create the response list 
response_ind_l = []

# the number of rolling 
till_end_bool = False
if till_end_bool == False:
    roll_num = 5 # if wanna go to the end specify in the loop

# the window size specification
train_size = 60 # days for training window
vali_size = 20 # days for validation window
test_size = 1 # days for testing window

# setup string
setup_str = '5_5_61_setI'

# option_type 
option_type = 'call'
#option_type = 'put'

for i in range(1, (1 + 2* (knots_K + degree_m + 1 - 1)+1)):  # num of alphas , last + 1 for range()
    response_ind_l.append('alpha' + str(i))

## response_ind 
#response_ind = '2_alpha1'

for response_ind_cur in response_ind_l[0:3]: # [0,6]
        
    #dat = pd.read_pickle('/work/LAS/cindyyu-lab/zc_gpu/gkx_all_inter_sic_v3')
    dat = pd.read_csv('/work/LAS/cindyyu-lab/zc_gpu/IV_real_surface_s_v2/surface_s_v2_alpha_dfs_' + option_type + '_'  + setup_str + '/'+ response_ind_cur +'.csv',index_col=0)
    if till_end_bool == True:  
        roll_num = dat.shape[0] - (train_size + vali_size + test_size) + 1 # if wanna go to the end specify in the loop
 
    test_date_l = dat['trade_date'][(train_size + vali_size):(train_size + vali_size + test_size * roll_num)].reset_index(drop = True)
    
    dat.reset_index(inplace = True, drop = True) # adjust index order
    dat.rename(columns={"response_value": "y"}, inplace = True)
    
    # drop the variables that are not predictors or response
    dat.drop(['trade_date','response_ind'], axis = 1,inplace = True)

   # dataframe to save results for current response_ind_cur
    temp_pred_save_XGB_df = pd.DataFrame()
     
    for x in range(1, roll_num + 1): # of repeated iterations, AT: start from 1
        print('response: ' + response_ind_cur)
        print('x value: ' + str(x))
    #    year_obs_count = pd.read_pickle('/work/LAS/cindyyu-lab/zc_gpu/GKX_real/year_obs_count.pkl') 
    #    init = period_gen(year_obs_count, x)[0]
    #    hori = period_gen(year_obs_count, x)[1] 
    #    tesp = period_gen(year_obs_count, x)[2] 
    #    starting_buffer = year_obs_count.iloc[:20 + x].sum() # fixed training window not expanding
        init = train_size # days for training window
        hori = vali_size # days for validation window
        tesp = test_size # days for testing window
        starting_buffer = 0 + (x-1)*tesp # fixed training window not expanding, everytime roll size of tesp forward
    
    ################################################# For real data
        dat_m1 = dat
        print(dat_m1.shape)
        # Rolling procedure: T = 180 (if taken as monthly data), N = 200 
        tscv = expanding_window(initial = init , horizon = hori, period = 100, test_p = tesp)
    
    #    tscv.split(dat_m1.iloc[:1102932,:]) # inefficient
    
        # Notcie that iloc does not include the endpoint ! .loc did 
    #    tscv.assign_ind(dat_m1, [list(range(init))], [list(range(init, init + hori))], [list(range(init + hori, init + hori + tesp))])
        tscv.assign_ind(dat_m1, [list(range(starting_buffer, init + starting_buffer))], [list(range(init + starting_buffer, init + hori + starting_buffer))], [list(range(init + hori + starting_buffer, init + hori + tesp + starting_buffer))])
        print('finish splitting')
        print(datetime.now())
    #    PCA
    #    result_pca.append(tscv.pca_tune(30,pp)) # max PCs
    
    # skopt 
    #    result_pca = tscv.pca_tune_skopt([1, 100], 100) # n_cores and max PCs 
    #    print('finish pca') 
    #    print(datetime.now())
    #    pd.DataFrame(result_pca).to_pickle('/work/LAS/cindyyu-lab/zc_gpu/GKX_real/result_GKX/pca_r/result_pca_single_' + str(x) + '.pkl')
    #    tscv.pred_save_pca().to_pickle('/work/LAS/cindyyu-lab/zc_gpu/GKX_real/result_GKX/pca_r/predictions_single_' + str(x) + '.pkl')
    
    #    PLS  
    #    result_pls.append(tscv.pls_tune(30,pp)) # n_cores and max PCs 
    
    # skopt 
    #    result_pls = tscv.pls_tune_skopt([1, 50], 100) # n_cores and max PCs 
    #    print('finish pls')
    #    print(datetime.now())   
    #    pd.DataFrame(result_pls).to_pickle('/work/LAS/cindyyu-lab/zc_gpu/GKX_real/result_GKX/pls_r/result_pls_single_' + str(x) + '.pkl')
    #    tscv.pred_save_pls().to_pickle('/work/LAS/cindyyu-lab/zc_gpu/GKX_real/result_GKX/pls_r/predictions_single_pls' + str(x) + '.pkl')
    
    
    #    Enet without Huber Loss 
    #    result_Enet.append(tscv.Enet_tune(list(10 ** (np.linspace(-2, 4, 61))), pp))
    
    #    result_Enet = tscv.Enet_tune_skopt([1e-1, 1], 100)
    #    print('finish Enet')
    #    print(datetime.now())
    #    pd.DataFrame(result_Enet).to_pickle('/work/LAS/cindyyu-lab/zc_gpu/GKX_real/result_GKX/result_Enet_longest.pkl')
    
    
    
    
    # Random forest skopt
        result_XGB = tscv.xgb_tune_skopt([2, 9], [1e-4, 1e-1], 100, 60) 
        print('finish XGB')
        print(datetime.now())
    
        temp_pred_df = tscv.pred_save_XGB()
        temp_pred_save_XGB_df = temp_pred_save_XGB_df.append(temp_pred_df)
#    pdb.set_trace()    
    # add trade_date column
    temp_pred_save_XGB_df.reset_index(drop = True, inplace = True)
    print(temp_pred_save_XGB_df)

    temp_pred_save_XGB_df.insert(0, 'trade_date', test_date_l)
    temp_pred_save_XGB_df.to_csv('/work/LAS/cindyyu-lab/zc_gpu/IV_real_surface_s_v2/result_remote_surface_s_v2_alphas_' + option_type + '_' + setup_str + '/pred_y/'+ response_ind_cur+ '_predictions_XGB'+ '.csv', index=False)

    end = time.time()
    print(end - start)
