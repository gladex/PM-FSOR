# PM-FSOR
Prognostic model based on a feature selection with orthogonal regression approach.

pca.R:
PCA analysis and visualization procedures for four LUAD GEO datasets.

FSOR_datapre.R:
Pretreat the gene expression data in TCGA and set initial parameters for FSOR.

FSOR_train.R:
It contains the main training process of FSOR. It may be time-consuming. Please split the experimental data small if necessary.

prognostic_model.R:
This involves further screening of genes based on cox regression results. C-index were considered in this procedure. 

lasso_cox.R:
The results were compared with THOSE of FSOR.

Contact:
Yuqi Wang (yukihhu@foxmail.com); Binhua Tang, PhD (bh.tang@outlook.com).
