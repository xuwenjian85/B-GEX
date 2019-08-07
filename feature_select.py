#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:05:37 2019

@author: xwj
"""
from multiprocessing import Pool
from functools import partial

import pandas as pd
import time, os
from sklearn import linear_model
from sklearn.feature_selection import SelectFromModel

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

