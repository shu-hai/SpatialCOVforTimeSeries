ptm<-proc.time()
library("MASS")   #for mvrnorm
library("flare")  #for CLIME
source('sugm_shu.R')

alpha=2
p=400#400

n=200

############# Approximate hyperbolic rate
#exact solve,optimal N is chosen by MAE
NN=1:9
error=numeric(length(NN))
k=log10(n)
x=1:n

f=x^(-alpha)
ind.const.posi=rep(1,length(NN)) #indicator for positive const

for(nn in 1:length(NN))
{
  N=NN[nn]
  
  beta=10^(k/N)
  
  
  A=matrix(0,N+1,N+1)
  temp=0:N
  for(i in temp)
  {
    A[(i+1),]=beta^(-temp*alpha)*exp(alpha)*exp(-alpha*beta^(i-temp))
  }
  Y=(beta^temp)^-alpha
  const=solve(A,Y)
  
  if(sum(const<0)>0)
  {
    ind.const.posi[nn]=-1
  }
  
  g=numeric(n)
  for(i in 0:N)
  {
    g=g+const[i+1]*beta^(-i*alpha)*exp(alpha)*exp(-alpha/beta^i*x)
  }
  
  
  
  error[nn]=sum(abs(f-g))
  
}

######
#optimal N
error=error*ind.const.posi
N=NN[which(error==min(error[error>=0]))]

beta=10^(k/N)
A=matrix(0,N+1,N+1)
temp=0:N
for(i in temp)
{
  A[(i+1),]=beta^(-temp*alpha)*exp(alpha)*exp(-alpha*beta^(i-temp))
}
Y=(beta^temp)^-alpha
const=solve(A,Y)


#####parameters transfer
w=numeric(N+1)
lambda=numeric(N+1)

for(i in 0:N)
{
  w[i+1]=const[i+1]*beta^(-i*alpha)*exp(alpha)
  lambda[i+1]=alpha/beta^i
}

a=sqrt(w*exp(-lambda))
rho=exp(-lambda)
#############
N=N+1 #includes index 0
###

##############  Model 1 simulated data
Omega=matrix(0,p,p) #precision matrix
Omega=0.6^abs(row(Omega)-col(Omega)) #AR(1) matrix with parameter 0.6

####generate Y and e vector
num=N*n
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) # seed for mvrnorm
mvn.vec=mvrnorm(n=num,mu=rep(0,p),Sigma=solve(Omega)) #num*p matrix

## X 
X=matrix(0,p,n)
## Y.1 
Y=mvn.vec[1:N,] # N*p matrix
## X.1
X[,1]=a%*%Y


t=N

for(k in 2:n)
{ 
  for(i in 1:N)
  {
    t=t+1
    Y[i,]=rho[i]*Y[i,]+sqrt(1-rho[i]^2)*mvn.vec[t,]
  }
  X[,k]=a%*%Y
}

rm(Y)
rm(mvn.vec)

########### CLIME estimator
##sample covariance matrix
S=cov(t(X))*(n-1)/n

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

##cross-validation
loss.tau=numeric(n.tau)
#H1=H2=10
n.blk=n/10
X_i=c((1:10)*n.blk,sample(n.blk:n,10))



 for(i in 1:20)
 {
  X_vld_e=X_i[i]
  X_vld_b=X_i[i]-(n.blk-1)
  
  X_vld=X[,X_vld_b:X_vld_e]# validation data
  
  X_trn_l=X_vld_b-(n.blk+1)
  X_trn_r=X_vld_e+(n.blk+1)
  
  if(X_trn_l<1){
    X_trn=X[,X_trn_r:n] 
  }
  else if(X_trn_r>n)
  {
    X_trn=X[,1:X_trn_l]
  }  
  else{
    X_trn=cbind(X[,1:X_trn_l],X[,X_trn_r:n])# training data
  }
  
  n.vld=dim(X_vld)[2]
  S.vld=cov(t(X_vld))*(n.vld-1)/n.vld
  
  n.trn=dim(X_trn)[2]
  S.trn=cov(t(X_trn))*(n.trn-1)/n.trn
  
  clime.trn=sugm_shu(S.trn, samplesize=n.trn, lambda=tau,method = "clime",max.ite = 1e3,perturb = TRUE,verbose = FALSE) # run clime
  
  
	for(j in 1:n.tau)
	{
		clime.icov.trn=clime.trn$icov[[j]]
		loss.tau[j]=loss.tau[j]+  sum(diag((  S.vld%*%clime.icov.trn - diag(1, p))^2))
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
clime=sugm_shu(S, samplesize=n, lambda=tau.opt,method = "clime",perturb = TRUE,verbose = FALSE) 
Omega.est=clime$icov[[1]]

delta=Omega-Omega.est

dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")

#check positive definiteness
#min.eigen=eigen(Omega.est)$values[p]
#pd=sum(min.eigen>0) #indicator of positive definiteness
#loss=0
##if(pd==1)
#{
#  loss=sum(diag(S%*%Omega.est))-log(det(Omega.est))
#}

time=proc.time()-ptm
output=c(dist.norm.2,dist.norm.F,dist.norm.1,time)#,min.eigen,pd,loss,time)
  

write.table(output,file=paste("model1_50_1e3_tr_0d01_n",n,"_p",p,"_alpha",alpha,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)
