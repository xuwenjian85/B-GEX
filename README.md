# B-GEX
Blood-based multi-tissue gene expression inference with Bayesian regression
## This is a quick demo
As a demo, we prepared a small expression TPM dataset of 2 target tissues and their associated blood samples. Each target tissue have 4 target genes. `input.zip` is demo data. The `demo.py` will extract blood features, build linear regression models and evaluate inference model performances.
```shell
unzip demo/input.zip -d demo/
python demo.py
```
## Preprocessing
Download two data files from [GTEx Portal](https://gtexportal.org/home):
* GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct
* GTEx_v7_Annotations_SampleAttributesDS.txt

Prepare feature file and target file. The preprocessing script prepare INPUT gene TPM files for each TISSUE. It will take ~ 40 mins CPU time and 20GB RAM.
```
python preprocess.py
```
## Feature selection
Make a series of feature sets with different alpha values. Train and evaluate baseline LSR model with these feature sets. This step will take several weeks of CPU time and better to do on HPC cluster. 
```
python feature_select.py
python baseline_model_train_valid.py
```
An optimal feature set which have the smallest 5-fold Cross-validation MAE can be determined by analyzing the results.  
## Model training and cross-validation with optimal feature set
The optimal feature set information is provided in 26tissues_samplesize.txt. This step will take a couple weeks of CPU time and need a HPC cluster. 
```shellscript
python 4_model_train_valid.py
```
