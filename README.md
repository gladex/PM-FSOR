# PM-FSOR
Construction of a prognostic model based on feature selection with orthogonal regression (FSOR).

>/Data/ contains:
* FSOR_result.Rdata: Weight of genes and iteration process of FSOR training.
* geodata.Rdata: Expression matrix and group list of four GEO datastes.
* hrgene_exp_1060.Rdata: Expression matrix of 1060 up-regulated genes in TCGA.
* unicox_results.Rdata: Univariate COX regression results of DEGs.

>/Results/ contains:
* BP_GO.csv: Results of biological process in GO function enrichment analysis.
* CC_GO.csv: Results of cellular component in GO function enrichment analysis.
* MF_GO.csv: Results of molecular function in GO function enrichment analysis.
* kegg.csv: Results of KEGG pathway enrichment analysis.
* selectgene_effectsize.csv: Combined effect size of up-regulated genes selected from four GEO datastes.
* selectgene_fsor_50_weight.csv: Weights of top 50 genes of FSOR.
* tcga_Up_DEM.csv: Detail information of up-regulated DEGs in TCGA.
* tcga_univariate COX analysis_log2(x+1).csv: Univariate COX regression results of DEGs.

>/Scripts/ contains: 
* pca.R: Implementation of PCA analysis and visualization procedures for four LUAD GEO datasets.
* FSOR_datapre.R: Preprocessing the gene expression data from TCGA and initializing parameters for FSOR.
* FSOR_train.R: It contains the main training process of FSOR. It may be time-consuming. Please split the experimental data small if necessary.
* prognostic_model.R: This involves further screening candidate genes based on Cox regression results. C-index were considered in this procedure. 
* lasso_cox.R: The results were compared with those based on the FSOR approach.

## Contact:
Yuqi Wang (yukihhu@foxmail.com); Binhua Tang, PhD (bh.tang@outlook.com).

## Citation:
Tang B., Wang Y., Chen Y., Li M., Tao Y. A Novel Early-Stage Lung Adenocarcinoma Prognostic Model Based on Feature Selection With Orthogonal Regression. Frontiers in Cell and Developmental Biology, 2021(8), [ArticleID:620746](https://www.frontiersin.org/articles/10.3389/fcell.2020.620746/full).
