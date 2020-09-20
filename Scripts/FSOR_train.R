rm(list = ls())
options(stringsAsFactors = F)
options(digits = 12)
setwd('./Data')
load('FSOR_training_data.Rdata')
library(Matrix)
i<-1
j<-1
k<-1
Theta<-Theta1
goal_fun<-NULL
goal_fun[1]<-opfun1[1]
while (opfun1[i]>0.5) {
  A<-Theta%*%X%*%H%*%t(X)%*%t(Theta)#1060*1060
  B<-Theta%*%X%*%H%*%t(Y)#1060*2
  A_ev<-eigen(A)
  alpha<-Re(A_ev$values[1])
  A_P<-alpha*Id-A
  W<-W1
  opfun1<-NULL
  opfun2<-NULL
  opfun1[1]<-sum(diag(t(t(W)%*%Theta%*%X%*%H-Y%*%H)%*%(t(W)%*%Theta%*%X%*%H-Y%*%H)))
  opfun2[1]<-opfun1[1]
  i<-1
  j<-1
  diff1<-opfun1[1]
  diff2<-opfun2[1]
  while (diff1[i]>0.001) {
    M<-2*A_P%*%W+2*B
    # S<-svd(M,nu=min(nrow(M),p[1]),nv=min(nrow(M),p[1]))
    S<-svd(M)
    W<-(S$u)%*%t(S$v)
    i<-i+1
    # opfun1[i]<-sum(diag(t(W)%*%A%*%W-2*t(W)%*%B))
    opfun1[i]<-sum(diag(t(t(W)%*%Theta%*%X%*%H-Y%*%H)%*%(t(W)%*%Theta%*%X%*%H-Y%*%H)))
    diff1[i]<-opfun1[i-1]-opfun1[i]
    cat('Update W:','Iteration=',i-1,"Diff=",diff1[i],"goal_fun=",opfun1[i],"\n")
  }  
  k<-k+1
  goal_fun[k]<-opfun1[i]
  Q<-(X%*%t(H)%*%t(X))*(W%*%t(W))
  s<-as.matrix(diag(2*X%*%H%*%t(Y)%*%t(W)))
  rho<-1.001
  lambda2<-0.00
  miu<-0.01
  theta<-diag(diag(x=rep(1/dim(X)[1],dim(X)[1])))
  v<-theta
  lambda1<-matrix((x=rep(0,dim(X)[1])),dim(X)[1],1)
  oned<-matrix((x=rep(1,dim(X)[1])),dim(X)[1],1)
  Theta<-diag(x=theta[1:dim(X)[1]],dim(X)[1])
  opfun2<-NULL
  opfun2[1]<-sum(diag(t(t(W)%*%Theta%*%X%*%H-Y%*%H)%*%(t(W)%*%Theta%*%X%*%H-Y%*%H)))
  diff2<-opfun2[1]
  j<-1
  while (diff2[j]>0.001) {
    E<-2*Q+miu*Id+miu*oned%*%t(oned)
    f<-miu*v+miu*oned-lambda2*oned-lambda1+s
    theta<-solve(E)%*%f
    # theta<-theta/sum(theta)
    mat<-theta+(1/miu)*lambda1
    v<-ifelse(mat<0,0,mat) #pos(theta+(1/miu)*lambda1)
    lambda1<-lambda1+miu*(theta-v)
    lambda2<-as.numeric(lambda2+miu*(t(theta)%*%oned-1))
    miu<-rho*miu
    Theta<-diag(x=theta[1:dim(X)[1]],dim(X)[1])
    j<-j+1
    opfun2[j]<-sum(diag(t(t(W)%*%Theta%*%X%*%H-Y%*%H)%*%(t(W)%*%Theta%*%X%*%H-Y%*%H)))
    diff2[j]<-opfun2[j-1]-opfun2[j]
    cat('Update theta:','Iteration=',j-1,"Diff=",diff2[j],"goal_fun=",opfun2[j],"\n")
  }
  if (k>15) {
    if ((goal_fun[k-5]-goal_fun[k]<0.5)) {
      break
    }
  }
}
save(v,goal_fun,file = 'FSOR_result.Rdata')