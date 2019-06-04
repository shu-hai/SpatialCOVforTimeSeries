ptm<-proc.time()
library("MASS")   #for mvrnorm
library("flare")  #for CLIME
source('sugm_shu.R')

p=400#400

n=200

##############  Model 1 i.i.d. simulated data
Omega=matrix(0,p,p) #precision matrix
Omega=0.6^abs(row(Omega)-col(Omega)) #AR(1) matrix with parameter 0.6

## X
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) # seed for mvrnorm
X=mvrnorm(n=n,mu=rep(0,p),Sigma=solve(Omega)) #n*p matrix

########### CLIME estimator
##sample covariance matrix
S=cov(X)*(n-1)/n

##candidate tau
n.tau=50
tau.min.ratio = 0.01# default in flare is 0.4
tau.max.tmp1 = min(max(S - diag(diag(S))), -min(S - diag(diag(S))))
tau.max.tmp2 = max(max(S - diag(diag(S))), -min(S - diag(diag(S))))
if(tau.max.tmp1 == 0){ 
  tau.max = tau.max.tmp2
}else{
  tau.max = tau.max.tmp1
} 
tau.min = tau.min.ratio * tau.max
tau = exp(seq(log(tau.max), log(tau.min), length = n.tau))

##Ordinary 10-fold cross-validation
loss.tau=numeric(n.tau)

n.blk=n/10
X_i=(1:10)*n.blk


  for(i in 1:10)
  {
    X_vld_e=X_i[i]
    X_vld_b=X_i[i]-(n.blk-1)
    
    X_vld=X[X_vld_b:X_vld_e,]# validation data
    
    X_trn_be=setdiff(1:200,X_vld_b:X_vld_e)
    
    X_trn=X[X_trn_be,]       # training data
    
    n.vld=dim(X_vld)[1]
    S.vld=cov(X_vld)*(n.vld-1)/n.vld
    
    n.trn=dim(X_trn)[1]
    S.trn=cov(X_trn)*(n.trn-1)/n.trn
    
	  clime.trn=sugm_shu(S.trn, samplesize=n.trn, lambda=tau,method = "clime",max.ite = 1e3,perturb = TRUE,verbose = FALSE) # run clime
	  
	  
	  for(j in 1:n.tau)
	  {
		  clime.icov.trn=clime.trn$icov[[j]]
		  loss.tau[j]=loss.tau[j]+ sum(diag((  S.vld%*%clime.icov.trn - diag(1, p))^2))
	  }
  }


rm(X)
rm(X_vld)
rm(X_trn)
rm(S.trn)
rm(S.vld)
rm(clime.trn)
rm(clime.icov.trn)

#optimal tau
loss.tau[which(loss.tau==-Inf)]=NaN
index.tau.opt=which(loss.tau==min(loss.tau,na.rm = TRUE))
tau.opt=tau[index.tau.opt][1]

#clime for the whole data
clime=sugm_shu(S, n.trn=n, lambda=tau.opt,method = "clime",perturb = TRUE,verbose = FALSE) 
Omega.est=clime$icov[[1]]

delta=Omega-Omega.est

dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")

#check positive definiteness
#min.eigen=eigen(Omega.est)$values[p]
#pd=sum(min.eigen>0) #indicator of positive definiteness
#loss=0
#if(pd==1)
#{
#  loss=sum(diag(S%*%Omega.est))-log(det(Omega.est))
#}
time=proc.time()-ptm
output=c(dist.norm.2,dist.norm.F,dist.norm.1,time)#,min.eigen,pd,loss,time)



write.table(output,file=paste("model1_50_1e3_tr_0d01_iid_n",n,"_p",p,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)


