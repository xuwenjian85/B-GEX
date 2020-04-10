#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 11:12:06 2019
# This is a demo
 
@author: xwj
"""
from multiprocessing import Pool
from functools import partial

import pandas as pd
import time, os
from sklearn import linear_model
from sklearn.feature_selection import SelectFromModel
import numpy as np

from sklearn.model_selection import GridSearchCV, cross_val_score, KFold

def f2(tgene, top_set, features, targets): # L1
    y = targets[tgene]
    i_range = {1:0.16,
               2:0.14,
               3:0.12,
               4:0.10, 
               5:0.08, 
               6:0.06, 
               7:0.04, 
               8:0.02}
    # different alpha(i) will different number of selected features.
    lmla = linear_model.Lasso(alpha=i_range[top_set], max_iter=10000).fit(features,y)
    lmla_idx = SelectFromModel(lmla, prefit=True).get_support(indices=True)
    return features.columns[lmla_idx]

#%%    
TISSUE_list = pd.read_table('demo/TISSUE_list_example.txt',index_col=0)
outdir = "demo/select_features/"
os.system("mkdir -p "+ outdir)
for TISSUE in TISSUE_list.SMTSD[TISSUE_list['BTsample'] > 100]:
    print(TISSUE)
    FEATURES = pd.read_table('demo/input/'+TISSUE+'_features.txt',index_col=0)
    TARGETS = pd.read_table('demo/input/'+TISSUE+'_targets.txt',index_col=0)
    num = TARGETS.shape[1]
    CORE = 2
    for TOP_SET in range(1,9):
        prefix = TISSUE + "_"+ str(TOP_SET) + "_" + str(num)
        print(prefix)
        print(time.ctime(),TISSUE, 'gene:', num, 'top:',TOP_SET,'cpu:',CORE, 'start')
        with Pool(CORE) as p:
            tgenes= TARGETS.columns[:num]
            res2 = p.map(partial(f2, top_set=TOP_SET, features = FEATURES, targets = TARGETS), tgenes)
        df = pd.DataFrame(res2).T
        df.columns = tgenes
        df.to_csv(outdir + prefix+".txt", sep= "\t", index=False)
        print(time.ctime(),TISSUE, 'gene:', num, 'top:',TOP_SET,'cpu:',CORE, 'finish')

#%%
# calc validation metrics per gene * NUM_TRIALS
def metric_worker(x, m, sf, f, t, pg, N, K):
    '''m=mymodel, sf=select_features, f=features, t=targets, pg=p_grid, N=NUM_TRIALS, K'''
    # get selected features, get x, y from data
    t_gene = t.columns[x]
    #print(time.ctime(), x, t_gene)
    X, y = f[sf[t_gene].dropna()], t[t_gene]
    mae_scores, mse_scores = np.zeros(N), np.zeros(N)
    # Loop for each trial
    for i in range(N):
        # Choose cross-validation techniques for the inner and outer loops, independently of the dataset.
        inner_cv = KFold(n_splits=K, shuffle=True, random_state=i)
        outer_cv = KFold(n_splits=K, shuffle=True, random_state=i)
        # Non_nested parameter search and scoring
        clf = GridSearchCV(m, param_grid=pg, cv=inner_cv, iid=True, scoring='neg_mean_absolute_error')
        # Nested CV with parameter optimization
        mae_scores[i] = \
        cross_val_score(clf, X, y, cv=outer_cv, scoring='neg_mean_absolute_error', pre_dispatch=1, n_jobs =1).mean()
        mse_scores[i] = \
        cross_val_score(clf, X, y, cv=outer_cv, scoring='neg_mean_squared_error', pre_dispatch=1, n_jobs =1).mean()
        #print('numtrials',i, 't_gene',t_gene,'return')
    return [- mae_scores.mean(), np.sqrt( - mse_scores.mean()) ]

#%%
# calc validation metrics && optimize hyper parameters for ONE learning model
def calc_model_metric(mymodel, select_features, features, targets, p_grid, num):
    NUM_TRIALS, K, CORE = 1, 5, 2
    tgenes= targets.columns[0:num]
    with Pool(CORE) as p:
        resList = p.map(partial(metric_worker, m=mymodel, sf=select_features, 
                            f=features, t=targets, pg=p_grid, N=NUM_TRIALS, K=K), range(num))
       
    metrics_per_gene = pd.DataFrame(resList).T
    metrics_per_gene.columns, metrics_per_gene.index = tgenes, ["MAE", "RMSE"]#['MAE','r2','RMSE']
    # multiprocessing END
        # calc model performance metrics - ( mean +/- std )
    vmean = np.mean(metrics_per_gene.loc['MAE',:][:num])
    vstd = np.std(metrics_per_gene.loc['MAE',:][:num],ddof=1)
    #MAE = str(round(vmean,DIGIT))+ u"\u00B1" + str(round(vstd,DIGIT))
    MAE = ("%.4f" % vmean) + u"\u00B1" + ("%.4f" % vstd) 
#    vmean = np.mean(metrics_per_gene.loc['r2',:][:num])
#    vstd= np.std(metrics_per_gene.loc['r2',:][:num],ddof=1)
#    r2 = ("%.4f" % vmean) + u"\u00B1" + ("%.4f" % vstd) 
    vmean = np.mean(metrics_per_gene.loc['RMSE',:][:num]) 
    vstd = np.std(metrics_per_gene.loc['RMSE',:][:num],ddof=1)
    RMSE = ("%.4f" % vmean) + u"\u00B1" + ("%.4f" % vstd) 
    # output metrics_per_gene into file
    metrics_per_gene.to_csv("metrics_per_gene.txt",sep= '\t', index=True)
    return [MAE, RMSE]

#%%
models, metrics = ['LSR','LSR-L1','LSR-L2','BayR'] , ['MAE','RMSE']
TopSetBest = pd.read_table("26tissues_samplesize.txt",index_col=0)
indir1 = "demo/input/"
indir2 = "demo/select_features/"
outdir = "demo/model_metric/"
os.system("mkdir -p "+ outdir)

#%%

for TISSUE in  TISSUE_list.SMTSD:
    features = pd.read_table(indir1 + TISSUE+"_features.txt",index_col=0)
    targets = pd.read_table(indir1 + TISSUE+"_targets.txt",index_col=0)
    num = targets.shape[1]
    
    #for TOP_SET in range(2,9): default 
    for TOP_SET in [TopSetBest.loc[TISSUE,'besttopset']]: # use this line to compute only the best topset 
        '''when run task need to resume from certain point ...., add this [if] statement'''
        prefix = TISSUE + "_"+ str(TOP_SET) + "_" + str(num)
        print(prefix)
        select_features = pd.read_table(indir2 + prefix+".txt", low_memory=False)
        key = 'LSR-L1'
        # table : final metric values 
        model_metric_table = pd.DataFrame(data = 0, 
                                          columns= [key + "_" + str(TOP_SET)],
                                          index= pd.MultiIndex.from_product([models, metrics]))
        
        # Ordinary least squares Linear Regression.
        print(time.ctime(), 'feature selection:',key,' + model: LSR')
        p_grid = {}
        mymodel = linear_model.LinearRegression()
        model_metric_table.loc[['LSR'],:] = \
            calc_model_metric(mymodel, select_features, features, targets, p_grid, num)
        os.rename("metrics_per_gene.txt", outdir + prefix+ "_LSR_metrics_per_gene.txt")
        
        # LASSO Linear Model trained with L1 prior as regularizerSO Linear Model trained with L1 prior as regularizer 
        print(time.ctime(), 'feature selection:',key,' + model: Lasso')
        p_grid = {'alpha': [1e-3, 0.01, 0.1, 1, 10]}
        mymodel = linear_model.Lasso(max_iter=100000)
        model_metric_table.loc[['LSR-L1'],:] = \
            calc_model_metric(mymodel, select_features, features, targets, p_grid, num)
        os.rename("metrics_per_gene.txt", outdir + prefix+ "_Lasso_metrics_per_gene.txt")

        # Linear least squares with l2 regularization.
        print(time.ctime(), 'feature selection:',key,' + model: Ridge')                
        mymodel = linear_model.Ridge()    
        p_grid =  {'alpha': [1e-3,0.01, 0.1, 1, 10,100,1000]}

        model_metric_table.loc[['LSR-L2'],:] = \
            calc_model_metric(mymodel, select_features, features, targets, p_grid, num)
        os.rename("metrics_per_gene.txt", outdir + prefix+ "_Ridge_metrics_per_gene.txt")
           
        # Bayesian Ridge Regression   
        print(time.ctime(), 'feature selection:',key,' + model: BayR')            
        mymodel = linear_model.BayesianRidge(n_iter=300)
        p_grid = {}
        model_metric_table.loc[['BayR'],:] = \
            calc_model_metric(mymodel, select_features, features, targets, p_grid, num)
        os.rename("metrics_per_gene.txt", outdir + prefix+ "_BayR_metrics_per_gene.txt")

        model_metric_table.to_csv(outdir + prefix+ ".txt",sep= '\t', index=True)
