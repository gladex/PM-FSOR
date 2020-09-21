# PM-FSOR
Construction of a prognostic model based on feature selection with orthogonal regression.




/Scripts/ contains: 
* pca.R: Implementation of PCA analysis and visualization procedures for four LUAD GEO datasets.
* FSOR_datapre.R: Preprocessing the gene expression data from TCGA and initializing parameters for FSOR.
* FSOR_train.R: It contains the main training process of FSOR. It may be time-consuming. Please split the experimental data small if necessary.
* prognostic_model.R: This involves further screening candidate genes based on Cox regression results. C-index were considered in this procedure. 
* lasso_cox.R: The results were compared with those based on the FSOR approach.

Contact:
Yuqi Wang (yukihhu@foxmail.com); Binhua Tang, PhD (bh.tang@outlook.com).
