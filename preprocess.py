#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 15:30:34 2019
@author: xwj
output feature selection txt with feature selection method LSR-L1
by iterate all tissues with > 100 available blood samples in GTex RNA seq TPM dataset
"""
import pandas as pd
import numpy as np
from sklearn import  preprocessing

#%%
gene_tpm = pd.read_table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct",skiprows=2)
samples = pd.read_table("GTEx_v7_Annotations_SampleAttributesDS.txt")
samples.insert(loc=0, column='SUBJID', value =np.repeat(0,samples.shape[0]))
'''add one column SUBJID'''
for i in range(samples.shape[0]):
    temp=samples.loc[i,'SAMPID'].split('-')
    samples.loc[i,'SUBJID']=temp[0]+'-'+temp[1]
samples = samples[samples.SMAFRZE == 'RNASEQ']

#%%
'''write blood and tissue RNAseq samples basic information'''
TISSUE_list = samples.loc[(samples.SMTSD != 'Whole Blood'), ['SMTS','SMTSD']]
TISSUE_list.drop_duplicates(inplace=True)

#%%
''' [take ~ 5 mins] ''' 
for t in TISSUE_list.index:
    TISSUE = TISSUE_list.SMTSD[t]
    target_id = samples[(samples.SMTSD == TISSUE ) & (samples.SMAFRZE == 'RNASEQ')].loc[:,['SUBJID','SAMPID']]
    feature_id = samples[(samples.SMTSD == 'Whole Blood') & (samples.SMAFRZE == 'RNASEQ')].loc[:,['SUBJID','SAMPID']]
    sample_overlap = list(set(target_id.SUBJID) & set(feature_id.SUBJID))
    feature_tpm = gene_tpm[["Name", "Description"] + list(feature_id[feature_id.SUBJID.isin(sample_overlap)].SAMPID)]
    target_tpm = gene_tpm[["Name", "Description"] + list(target_id[target_id.SUBJID.isin(sample_overlap)].SAMPID)]
    # fill tissue sample number statistics
    TISSUE_list.loc[t,'Bsample'] = feature_id.shape[0]
    TISSUE_list.loc[t,'Tsample'] = target_id.shape[0]
    TISSUE_list.loc[t,'BTsample'] = len(sample_overlap)
    TISSUE_list.loc[t,'Tgenes'] = target_tpm.shape[0]
    TISSUE_list.loc[t,'TgenesClean'] = target_tpm.loc[target_tpm.all(axis=1)].shape[0]
    TISSUE_list.loc[t,'Bgenes'] = feature_tpm.shape[0]
    TISSUE_list.loc[t,'BgenesClean'] = feature_tpm.loc[feature_tpm.all(axis=1)].shape[0]

#%%    
TISSUE_list = TISSUE_list[(TISSUE_list.BTsample > 100) & 
                          (TISSUE_list.SMTSD.values != 'Cells - Transformed fibroblasts') &
                          (TISSUE_list.SMTSD.values != 'Cells - EBV-transformed lymphocytes')]
TISSUE_list.index = range(TISSUE_list.shape[0])
TISSUE_list.to_csv('TISSUE_list.txt', sep='\t')

#%%    
'''
prepare INPUT gene TPM files for the TISSUEs [take ~ 30 mins] 
'''
for TISSUE in TISSUE_list.SMTSD:
    print(TISSUE) 
    target_id = samples[(samples.SMTSD == TISSUE) & (samples.SMAFRZE == 'RNASEQ')].loc[:,['SUBJID','SAMPID']]
    feature_id = samples[(samples.SMTSD == 'Whole Blood') & (samples.SMAFRZE == 'RNASEQ')].loc[:,['SUBJID','SAMPID']]
    sample_overlap = list(set(target_id.SUBJID) & set(feature_id.SUBJID))
    feature_tpm = gene_tpm[["Name", "Description"] + list(feature_id[feature_id.SUBJID.isin(sample_overlap)].SAMPID)]
    target_tpm = gene_tpm[["Name", "Description"] + list(target_id[target_id.SUBJID.isin(sample_overlap)].SAMPID)]
    '''data cleaning: remove incomplete columns of genes (which are not expressed in all Samples)'''
    feature_tpm = feature_tpm.loc[feature_tpm.all(axis=1)]
    target_tpm = target_tpm.loc[target_tpm.all(axis=1)]
    '''---------------------!! checkpoint 1----------------------------------------'''
    if feature_tpm.all().all() & (not feature_tpm.empty) & target_tpm.all().all() & (not target_tpm.empty):
        pass
    else:
        print(TISSUE, ':clean TPM error')
    '''-----------------------------------------------------------------------------'''
    '''reshaping: [Gene x Sample] -transpose-> [Sample x Gene]'''
    target_tpm_T = pd.DataFrame(data=0,index=target_tpm.columns[2:], columns=target_tpm.Name)
    gene_total = target_tpm.shape[0]
    step,i,j= 1000,0,0
    while i < gene_total:
        j = i + step if i + step < gene_total else gene_total
        target_tpm_T.iloc[:,i:j] = target_tpm.iloc[i:j,2:].transpose().values
        i=j
    
    feature_tpm_T = pd.DataFrame(data=0,index=feature_tpm.columns[2:], columns=feature_tpm.Name)
    gene_total = feature_tpm.shape[0]
    step,i,j= 1000,0,0
    while i < gene_total:
        j = i + step if i + step < gene_total else gene_total
        feature_tpm_T.iloc[:,i:j] = feature_tpm.iloc[i:j,2:].transpose().values
        i=j
    '''-----------------------------------------------------------------------------'''
    '''TPM value scaling：Standardize Per Gene mean=0，variance=1 '''
    data = preprocessing.StandardScaler().fit_transform(feature_tpm_T)
    feature_tpm_T_scaled = pd.DataFrame(data,index=feature_tpm_T.index,columns=feature_tpm_T.columns)
    
    data = preprocessing.StandardScaler().fit_transform(target_tpm_T)
    target_tpm_T_scaled = pd.DataFrame(data,index=target_tpm_T.index,columns=target_tpm_T.columns)
    del data, i, j, step,target_id, feature_id, sample_overlap, gene_total, feature_tpm, target_tpm
    '''-----------------------------------------------------------------------------'''
    '''output clean data as *.txt'''
    features, targets = feature_tpm_T_scaled, target_tpm_T_scaled
    #targets, features = feature_tpm_T_scaled, target_tpm_T_scaled
    features.to_csv('input/'+TISSUE+'_features.txt',sep='\t')
    targets.to_csv('input/'+TISSUE+'_targets.txt', sep='\t')
