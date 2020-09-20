rm(list = ls())
options(stringsAsFactors = F)
options(digits = 12)
setwd('./Data')
library(rstiefel)
library(Matrix)
library(survival)
library(magrittr)
library(tidyverse)
library(stringr)
load('gene_exp_1060.Rdata')
#load data
X<-hrgene_exp
X<-log2(X+1)
# X<-scale(X,center=apply(X,2,median),scale=apply(X,2,function(x){sum(abs(x-mean(x)))/(length(x))}))
Y<-LUAD_clinicaldata[,c("time_year","event")]
Y<-t(Y)
Y<-matrix(Y,nrow = 2,ncol = 479)
####calculate coefficient
Theta1<-diag(x=rep(1/dim(X)[1],dim(X)[1]))#1060*1060
Theta<-Theta1
In<-diag(x=rep(1,dim(Y)[2]))#
Id<-diag(x=rep(1,dim(X)[1]))#1060*1060
In_1<-matrix(rep(1/dim(Y)[2],479),nrow = 479,ncol = 479)
H<-In-In_1
A<-Theta%*%X%*%H%*%t(X)%*%t(Theta)#1060*1060
B<-Theta%*%X%*%H%*%t(Y)#1060*2
A_ev<-eigen(A)
alpha<-Re(A_ev$values[1])
set.seed(88)
W1<-matrix(rnorm(dim(X)[1]*2),dim(X)[1],2)
W1<-rmf.matrix(W1)
W<-W1
A_P<-alpha*Id-A
opfun1<-NULL
opfun2<-NULL
opfun1[1]<-sum(diag(t(t(W)%*%Theta%*%X%*%H-Y%*%H)%*%(t(W)%*%Theta%*%X%*%H-Y%*%H)))
opfun2[1]<-opfun1[1]
i<-1
j<-1
diff1<-opfun1[1]
diff2<-opfun2[1]
save(X,Y,W1,In,Id,H,Theta1,opfun1,opfun2,file = 'FSOR_training_data.Rdata')

