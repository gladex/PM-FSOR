rm(list = ls())
options(stringsAsFactors = F)
options(digits = 12)
setwd('./Data')
library(rstiefel)
library(Matrix)
library(survival)
library(ggplot2)
library(magrittr)
library(tidyverse)
library(stringr)
load('FSOR_result.Rdata')
#top50 genes
selectgene<-v[order(v,decreasing = TRUE)[1:50],]
load('unicox_results.Rdata')
# results<-results[results$`P-value`<0.05,]
results<-results[results$beta>0,]
results<-results[results$`P-value`<0.2,]
selectgenenames<-data.frame(names(selectgene))
names(selectgenenames)<-c("ID")
selectgenenames<-selectgenenames[selectgenenames$ID%in%rownames(results),]
selectgenenames<-data.frame(selectgenenames)
names(selectgenenames)<-c("ID")
Expr<-t(Expr)
Expr<-Expr[rownames(Up_DEM_expr),]
# selectgene_nameid<-selectgenenames%>%
#   separate(ID,into = c('GeneSymbol','GeneID'),sep = "[|]",remove = FALSE)
# write.csv(selectgene_nameid,file = 'selectgene_fsor_50_1.csv',fileEncoding = 'utf8')

# RFcolnames %>% 
#   separate(ID, into = c("Gene2", "IDnum2"),sep = "[.]",remove = FALSE) %>% 
#   head()
# selectgenenames2<-selectgenenames %>%
#   separate(ID, into = c("Gene2", "IDnum2"),sep = "[.]",remove = FALSE)
# selectgenenames2<-selectgenenames2 %>%
#   unite(ID_new, Gene2:IDnum2, sep = "|")
selectgene_expr<-Expr[,selectgenenames$ID]
selectgene_expr<-log2(selectgene_expr+1)
selectgene_expr<-data.frame(selectgene_expr)
mySurv<-with(LUAD_clinicaldata,Surv(time_year,event))
###ppi results
string_32<-read.csv('./Results/string_interactions_fsor_50_1.tsv default node.csv')
string_32gene<-selectgene_nameid[selectgene_nameid$GeneSymbol%in%string_32$name,]
geneweight_32<-geneweight_50[geneweight_50$GeneSymbol%in%string_32$name,]
string_32gene<-string_32gene %>%
     unite(symbolid, GeneSymbol:GeneID,sep = ".")
#######multicox regression#######
selectgenenames3_expr<-selectgene_expr[,string_32gene$symbolid]
muli_COX4<-coxph(mySurv ~ ., data=selectgenenames3_expr) 
summary(muli_COX4)
ph_hypo_multi<-cox.zph((muli_COX4))
ph_hypo_table<-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
formula_for_multivariate<-as.formula(paste0('mySurv~',paste(rownames(ph_hypo_table[ph_hypo_table[,3]>0.05,]),sep='',collapse = '+')))
muli_COX5<-coxph(formula_for_multivariate,data=selectgenenames3_expr)
summary(muli_COX5)
# step.multi_COX=step(muli_COX5,direction = "both")
# step.multi_COX 
beta<-muli_COX5[1]$coefficients
beta<-cbind(names(beta),beta)
beta<-data.frame(beta)
names(beta)<-c('gene','coefficients')
beta_gene<-beta$gene
# beta_gene<-beta[beta$coefficients>-0.12,]$gene
# beta_gene<-names(step.multi_COX$coefficients[step.multi_COX$coefficients>0])
# beta_gene<-names(step.multi_COX$coefficients)
# formula_multi_COX6<-as.formula(paste0('mySurv~',paste(beta_gene,sep='',collapse = '+')))
# multi_COX6<-coxph(formula_multi_COX6,data = selectgenenames3_expr)
# summary(multi_COX6)
# step.multi_COX=step(multi_COX6,direction = "both")
# step.multi_COX
# beta_gene<-names(step.multi_COX$coefficients)
#####stepwise c-index#######
# beta_gene<-beta[beta$coefficients>0,]$gene

#####select genes via C-index#######
c_result<-NULL
bestresult<-NULL
k<-1
formula_multi_COX6<-as.formula(paste0('mySurv~',paste(beta_gene,sep='',collapse = '+')))
multi_cox<-coxph(formula_multi_COX6,data = selectgenenames3_expr)
bestresult[1]<-multi_cox$concordance[6]
while(length(beta_gene)>0){
  formula_mcox<-as.formula(paste0('mySurv~',paste(beta_gene,sep='',collapse = '+')))
  multi_cox<-coxph(formula_mcox,data = selectgenenames3_expr)
  print(paste0('C-index:',multi_cox$concordance[6]),quote=F,row.names=F)
  print(multi_cox$formula,quote=F,row.names=F)
  j<-length(beta_gene)
  k<-k+1
  c_result<-NULL
  for (i in 1:j) {
  beta_gene_new<-beta_gene[-i]
  formula_mcox<-as.formula(paste0('mySurv~',paste(beta_gene_new,sep='',collapse = '+')))
  multi_cox<-coxph(formula_mcox,data = selectgenenames3_expr)
  c_result[i]<-multi_cox$concordance[6]
  }
  result_gene_c<-cbind(beta_gene,c_result)
  result_gene_c_order<-result_gene_c[order(result_gene_c[,2],decreasing = T),]
  beta_gene_order<-matrix(result_gene_c_order[,1])
  beta_gene_m<-apply(beta_gene_order,1,function(x){paste0("-",x)})
  result_gene_c_order_p<-data.frame(cbind(beta_gene_m,result_gene_c_order[,2]))
  names(result_gene_c_order_p)<-c('gene','C-index')
  result_gene_c_order_p$`C-index`<-as.numeric(result_gene_c_order_p$`C-index`)
  bestresult[k]<-result_gene_c_order_p$`C-index`[1]
  print(result_gene_c_order_p,quote = F,row.names=F)
  if(bestresult[k]-bestresult[k-1]<0){
    break
  }
  beta_gene<-beta_gene_order[-1]
}
formula_mcox<-as.formula(paste0('mySurv~',paste(beta_gene,sep='',collapse = '+')))
step.multi_COX<-coxph(formula_mcox,data = selectgenenames3_expr)
summary(step.multi_COX)
#########risk group###########
# step.multi_COX=step(step.multi_COX ,direction = "both")
# step.multi_COX #查看结果,最后筛选
# genenames_8<-names(step.multi_COX$coefficients)
# save(genenames_8,file = 'genenames_8.Rdata')
RiskScore<-predict(step.multi_COX,type = "risk",newdata =selectgenenames3_expr)
RiskScore 
risk_group<-ifelse(RiskScore>=median(RiskScore),'high','low')
risk_group
######
library(stringr)
str(LUAD_clinicaldata)
dim(LUAD_clinicaldata)
class(names(RiskScore))
all(substr(names(RiskScore),1,12)==substr(rownames(LUAD_clinicaldata),1,12))
all(substr(names(risk_group),1,12)==substr(rownames(LUAD_clinicaldata),1,12))
KM.input<-cbind(LUAD_clinicaldata[,c("event","time_year")],RiskScore,risk_group)
library(survival)
library(survminer)
str(KM.input)
fit<-survfit(Surv(time_year,event) ~ risk_group, data=KM.input)
print(fit)
summary(fit) 
summary(fit)$table 
#######
tiff(file="KM.tiff", width = 13, height =13, units ="cm",
     compression="lzw", bg="white", res=600)
ggsurvplot(fit, 
           size=1.2, pval=T, conf.int=F, 
           conf.int.style="step",
           #fun = "cumhaz",
           risk.table= T, #"abs_pct", # T 
           risk.table.fontsize=4,
           risk.table.y.text.col= T,
           risk.table.y.text= F,
           risk.table.col = "strata",
           risk.table.height = 0.22,
           tables.theme = theme_cleantable(),
           
           legend=c(0.8,0.8), legend.title="Strata",
           xlab= "Time (years)", ylab="Survival Probability",
           surv.median.line= "hv", 
           ncensor.plot= F,
           
           legend.labs= c("High risk","Low risk"),
           font.main=c(12,"bold","black"),
           font.x=c(12,"bold","black"),
           font.y=c(12,"bold","black"),
           font.tickslab=c(8,"bold","black"),
           palette=c("#E31A1C","#1D91C0"),
           #palette = "Dark2",
           #palette = c("#E7B800", "#2E9FDF"), 
           ggtheme= theme_bw())
#ggsave(file=paste("Fig.X.LUAD",".pdf", sep=""))
# width=4.14,height=6.14
dev.off()

#######ROC curve########
library(timeROC)
time_ROC_input<-KM.input
time_ROC<-timeROC(T=time_ROC_input$time_year,
                  delta=time_ROC_input$event, 
                  marker=time_ROC_input$RiskScore, 
                  cause=1, 
                  weighting = "marginal", 
                  times = c(1,2,3,4),
                  ROC=TRUE,
                  iid=TRUE
)
time_ROC 
summary(time_ROC) 
time_ROC$AUC
library(ggplot2)
time_ROC$TP
summary(time_ROC)
time_ROC.res<-data.frame(
                        TP_1year=time_ROC$TP[,1], 
                        FP_1year=time_ROC$FP[,1],
                        TP_2year=time_ROC$TP[,2], 
                         FP_2year=time_ROC$FP[,2],  
                         TP_3year=time_ROC$TP[,3],  
                         FP_3year=time_ROC$FP[,3], 
                         TP_4year=time_ROC$TP[,4],
                         FP_4year=time_ROC$FP[,4]) 
time_ROC$AUC
TimeROC_plot<-ggplot()+
  geom_line(data=time_ROC.res,aes(x=FP_1year,y=TP_1year),size=1,color="green")+
  geom_line(data=time_ROC.res,aes(x=FP_2year,y=TP_2year),size=1,color="red")+
  geom_line(data=time_ROC.res,aes(x=FP_3year,y=TP_3year),size=1,color="blue")+
  geom_line(data=time_ROC.res,aes(x=FP_4year,y=TP_4year),size=1,color="black")+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",size=1, linetype = 2 )+
  ggtitle('ROC curve of FSOR-COX')+
  theme_bw()+
  annotate("text",x=0.75,y=0.35,size=4.5,
           label=paste0("AUC at 1 years = ",round(time_ROC$AUC[[1]],3)),color="green")+
  annotate("text",x=0.75,y=0.25,size=4.5,
           label=paste0("AUC at 2 years = ",round(time_ROC$AUC[[2]],3)),color="red")+
  annotate("text",x=0.75,y=0.15,size=4.5,
           label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[3]],3)),color="blue")+
  annotate("text",x=0.75,y=0.05,size=4.5,
           label=paste0("AUC at 4 years = ",round(time_ROC$AUC[[4]],3)),color="black")+
  labs(x="False positive rate",y="True positive rate")+
  theme(plot.title  = element_text(hjust = 0.5,face = 'bold'),
    axis.text=element_text(face="bold", size=11,  color="black"),
        axis.title=element_text(face="bold", size=14, color="black")) 
p1<-TimeROC_plot 


