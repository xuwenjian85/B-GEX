# B-GEX
Blood-based multi-tissue gene expression inference with Bayesian regression.

The inference models are tissue-specific, so we built an independent dataset for each tissue. 

## This is a quick demo(to do)
As a demo, we prepared a small expression TPM dataset of 2 target tissues and their associated blood samples. Each target tissue have 4 target genes. `input.zip` is demo data. The `demo.py` will extract blood features, build linear regression models and evaluate inference model performances.
```shell
unzip demo/input.zip -d demo/
python demo.py
```
## Preprocessing
Download data files from [GTEx Portal](https://gtexportal.org/home) and gencode:
* GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct
* GTEx_v7_Annotations_SampleAttributesDS.txt
* gencode.v30lift37.annotation.gtf

Collect pseudogenes in a file. Then remove pseudogene from gene_tpm
```shell
awk -F "\t" '$3 == "gene" && $9~/pseudogene/ {print}' gencode.v30lift37.annotation.gtf | cut -f9 | cut -d" " -f2,4 | sed 's/[";]//g' | sed 's/[\. ]/\t/g' | cut -f1,3 > gencode.v30lift37.pseudogene.txt

cut -f1 gencode.v30lift37.pseudogene.txt | grep -v -f - GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct > GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.no_pseudogene
```
## Main script

The rest of main work can be done in B-GEX_v2.ipynb on a linux x86_64 server (2 x Intel(R) Xeon(R) Gold 6144 CPU 3.50GHz, 128 Gb RAM). 

## How it works?

### Feature selection

Feature is the blood  gene expression profile. Target is the tissue gene expression profile.

First, we compute coefficient of variation of each features and rank them in descending order. We select those top 10% features to construct a low dimensional new subsets.  

Second, we assume that one feature gene is important when the absolute of cosine similarity between vector. We use cosine similarity to reduce to [5,10,15,20,...] features per target gene. 

Train and evaluate baseline LSR model with these feature sets. We found 10 features is an optimal feature set for most tissues. 

### Model training with optimal feature set
The optimal feature set used to make a linear regression model for each gene in each tissue.

The models’ performance is further analyzed at gene level with MAE, r(Pearsons correlation coefficient), RMSE, etc.

For convenience, we separate genes into predictable and non-predictable genes according to empirical quality check results.          

> Rule 1, if target gene’s MAE < 0.7 and r > 0.3, the gene is predictable. 
>
> Rule 2. If the ratio of predictable genes to total genes > 0.2, the tissue is predictable.   

B-GEX outputs inferred gene expression values from blood gene expression profiles together with empirical quality check marks (MAE, RMSE, r, and predictability). Users can filter the output by these marks and should be cautious when the genes of interest fall into the non-predictable category.



