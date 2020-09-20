rm(list = ls()) 
options(stringsAsFactors = F)
setwd('./Data')
library(GEOquery)
library(stringr)
library(ggplot2)
library(ggthemes)
library(grid)
########data prepare######
#load geodata
load('geodata.Rdata')
dat1=exprs(gset1) 
dim(dat1)
dat2=exprs(gset2)
dat3=exprs(gset3)
dat4=exprs(gset4)
dat1[1:4,1:4]
dat1=log2(dat1+1)
dat2=log2(dat2+1)
dat3=log2(dat3+1)
dat4=log2(dat4+1)
# boxplot(dat,las=2)
#pdata
pd1=pData(gset1) 
pd2=pData(gset2)
pd3=pData(gset3)
pd4=pData(gset4)
#grouplist
group_list1<-ifelse(str_detect(pd1$characteristics_ch1.3,"histology: Adenocarcinoma"),"LUAD",
            ifelse(str_detect(pd1$characteristics_ch1.3,"histology: --"),"Normal","Others"))
table(group_list1)
group_list2<-ifelse(str_detect(pd2$source_name_ch1,"Lung adenocarcinoma"),"LUAD","Normal")
table(group_list2)
group_list3<-ifelse(str_detect(pd3$characteristics_ch1.3,"histology: adeno"),"LUAD",
                    ifelse(str_detect(pd3$characteristics_ch1.3,"histology: normal lung"),"Normal","Others"))
table(group_list3)
group_list4<-ifelse(str_detect(pd4$characteristics_ch1,"histology: Non-malignant"),"Normal","LUAD")
table(group_list4)
#######pca analysis########
#dat1
dat1=t(dat1)
dat1=as.data.frame(dat1)
dat1=cbind(dat1,group_list1) 
dat1=subset(dat1,group_list1=="LUAD"|group_list1=="Normal")
dfpca1<-prcomp(dat1[1:ncol(dat1)-1])
dfpcs1<-data.frame(dfpca1$x,Species=dat1[ncol(dat1)])
head(dfpcs1,3)
#dat2
dat2=t(dat2)
dat2=as.data.frame(dat2)
dat2=cbind(dat2,group_list2)
dfpca2<-prcomp(dat2[1:ncol(dat2)-1])
dfpcs2<-data.frame(dfpca2$x,species=dat2[ncol(dat2)])
#dat3
dat3=t(dat3)
dat3=as.data.frame(dat3)
dat3=cbind(dat3,group_list3) 
dat3=subset(dat3,group_list3=="LUAD"|group_list3=="Normal")
dfpca3<-prcomp(dat3[1:ncol(dat3)-1])
dfpcs3<-data.frame(dfpca3$x,Species=dat3[ncol(dat3)])
#dat4
dat4=t(dat4)
dat4=as.data.frame(dat4)
dat4=cbind(dat4,group_list4)
dfpca4<-prcomp(dat4[1:ncol(dat4)-1])
dfpcs4<-data.frame(dfpca4$x,species=dat4[ncol(dat4)])
####plot(dfpca1$x[,1], dfpca$x[,2])
#######plot###########
dfpcs1p=dfpcs1[order(dfpcs1[,"group_list1"]),]
dfpcs2p=dfpcs2[order(dfpcs2[,"group_list2"]),]
dfpcs3p=dfpcs3[order(dfpcs3[,"group_list3"]),]
dfpcs4p=dfpcs4[order(dfpcs4[,"group_list4"]),]
# ggplot(data=NULL,aes(dfpcs1$PC1,dfpcs1$PC2))+geom_point(data=dfpcs1,aes(colour=dfpcs1$group_list1))
# +geom_point(data=dfpcs2,aes(x=dfpcs1$PC1,y=dfpcs1$PC2))
p<-ggplot(dfpcs1p,aes(x=PC1,y=PC2))+
  geom_point(mapping=aes(shape="GSE32036",colour=dfpcs1p$group_list1),size=4)+
  geom_point(data=dfpcs2p,mapping=aes(x=PC2,y=PC1,shape="GSE32867",colour=dfpcs2p$group_list2),size=4)+
  geom_point(data=dfpcs3p,mapping=aes(shape="GSE33532",colour=dfpcs3p$group_list3),size=4)+
  geom_point(data=dfpcs4p,mapping=aes(x=PC2,y=PC1,shape="GSE75037",colour=dfpcs4p$group_list4),size=4)+
  labs(title="PCA for combined samples",shape="GEO dataset",colour='Type')+
  theme_bw()
p
ggsave(p,"pca.jpg",dpi=600)
# ggplot(dfpcs2p,aes(x=PC1,y=PC2))+geom_point(mapping=aes(colour=group_list2))+theme_bw()
# ggplot(dfpcs3p,aes(x=PC1,y=PC2))+geom_point(mapping=aes(colour=group_list3))+theme_bw()
# ggplot(dfpcs4p,aes(x=PC1,y=PC2))+geom_point(mapping=aes(colour=group_list4))+theme_bw()
# ggplot(data=NULL,aes(dfpcs1$PC1,dfpcs1$PC2))+theme_bw()
# plot(dfpca1$x[,1], dfpca1$x[,2],col="red")
# par(new=TRUE)
# plot(dfpca2$x[,1], dfpca2$x[,2])
# library("FactoMineR")
# library("factoextra") 
# # The variable group_list (index = 54676) is removed
# # before PCA analysis
# dat.pca <- PCA(dat1[,-ncol(dat1)], graph = FALSE)#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
# fviz_pca_ind(dat.pca,
#              geom.ind = "point", # show points only (nbut not "text")
#              col.ind = dat1$group_list, # color by groups
#              # palette = c("#00AFBB", "#E7B800"),
#              addEllipses = TRUE, # Concentration ellipses
#              legend.title = "Groups"
# )
# fviz_pca_ind(dat.pca,
#                  geom.ind = "text", # show points only (nbut not "text")
#                  col.ind = dat1$group_list, # color by groups
#                  # palette = c("#00AFBB", "#E7B800"),
#                  addEllipses = TRUE, # Concentration ellipses
#                  legend.title = "Groups"
# )
