rm(list = ls())
options(stringsAsFactors = F)
setwd('./Data')
library('ggplot2')
library('glmnet')
library('survival')
set.seed(123)
load('unicox_results.Rdata')
load('hrgene_exp_1060.Rdata')
hrgene_exp<-t(hrgene_exp)
hrgene_exp[1:6,1:6]
hrgene_exp<-log2(hrgene_exp+1)
hrgene_exp[1:6,1:6]
mySurv<-with(LUAD_clinicaldata,Surv(time_year,event))
hrgene_exp<-hrgene_exp[rownames(Up_DEM_expr),]
y<-LUAD_clinicaldata[,c("time_year","event")]
names(y)<-c('time','status')
y$time<-as.double(y$time)
y$status<-as.double(y$status)
y<-as.matrix(survival::Surv(y$time,y$status))
fit<-glmnet(hrgene_exp,y,family = 'cox')
plot(fit)
lasso_fit<-cv.glmnet(hrgene_exp,y,family='cox',type.measure = 'deviance')
plot(lasso_fit)
coefficient<-coef(lasso_fit,s=lasso_fit$lambda.min)
# coefficient<-coef(lasso_fit,s=lasso_fit$lambda.1se)
Active.Index<-which(as.numeric(coefficient)!=0)
active.coefficients<-as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox<-rownames(coefficient)[Active.Index]
lassogeneExpr<-hrgene_exp[,sig_gene_multi_cox]
save(lassogeneExpr,file = 'lassogene.Rdata')
########
mySurv<-with(LUAD_clinicaldata,Surv(time_year,event))
lassogeneExpr<-data.frame(lassogeneExpr)
multi_COX<-coxph(mySurv ~ ., data=lassogeneExpr) 
summary(multi_COX)
ph_hypo_multi<-cox.zph((multi_COX))
ph_hypo_table<-ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
formula_for_multivariate<-as.formula(paste0('mySurv~',paste(rownames(ph_hypo_table[ph_hypo_table[,3]>0.05,]),sep='',collapse = '+')))
muli_COX2<-coxph(formula_for_multivariate,data=lassogeneExpr)
summary(muli_COX2)
# step.multi_COX=step(muli_COX2,direction = "both")
step.multi_COX<-muli_COX2
####
RiskScore<-predict(step.multi_COX,type = "risk",newdata =lassogeneExpr)
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
#
KMsurvival_plot<-ggsurvplot(fit,pval = TRUE, 
                            conf.int = TRUE,
                            conf.int.style = "step", 
                            risk.table = "abs_pct",  
                            risk.table.y.text.col = T,
                            risk.table.y.text = FALSE,
                            xlab = "Time in years",  
                            surv.median.line = "hv",
                            ncensor.plot = FALSE, 
                            legend.labs =
                              c("high risk", "low risk"),    
                            palette = c("#3030FF", "#FF3939"), 
                            ggtheme = theme_light()
)
KMsurvival_plot
####ROC#####
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
  ggtitle('ROC curve of LASSO-COX')+
  theme_bw()+
  annotate("text",x=0.75,y=0.35,size=8,
           label=paste0("AUC at 1 years = ",round(time_ROC$AUC[[1]],3)),color="green")+
  annotate("text",x=0.75,y=0.25,size=8,
           label=paste0("AUC at 2 years = ",round(time_ROC$AUC[[2]],3)),color="red")+
  annotate("text",x=0.75,y=0.15,size=8,
           label=paste0("AUC at 3 years = ",round(time_ROC$AUC[[3]],3)),color="blue")+
  annotate("text",x=0.75,y=0.05,size=8,
           label=paste0("AUC at 4 years = ",round(time_ROC$AUC[[4]],3)),color="black")+
  labs(x="False positive rate",y="True positive rate")+
  theme(plot.title  = element_text(hjust = 0.5,face = 'bold',size = 20),
        axis.text=element_text(face="bold", size=11,  color="black"),
        axis.title=element_text(face="bold", size=15, color="black")) 
p2<-TimeROC_plot 
p2
