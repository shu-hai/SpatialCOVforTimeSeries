ptm<-proc.time()
library("MASS")   #for mvrnorm
library("QUIC")  #for spice


alpha=2
p=200#400

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

############# Model band simulated data
####construct Omega matrix
Omega=matrix(0,p,p)
for(i in 1:(p-2))
{
    
    Omega[i,i]=0.5
    Omega[i,(i+1)]=0.6
    Omega[i,(i+2)]=0.3
    
}
Omega[p-1,p-1]=0.5
Omega[p,p]=0.5
Omega[p-1,p]=0.6
Omega=Omega+t(Omega)



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

########### spice estimator
##sample covariance matrix
S=cov(t(X))*(n-1)/n

W_inv=diag(diag(S)^(-0.5))

R=W_inv%*%S%*%W_inv 
##candidate lambda, using huge package default lambda
n.lambda=50
lambda.min.ratio = 0.01# default in huge is 0.1
lambda.max = max(max(R - diag(p)), -min(R - diag(p)))
lambda.min = lambda.min.ratio * lambda.max
lambda = exp(seq(log(lambda.max), log(lambda.min), length = n.lambda))

##cross-validation
loss.lambda=numeric(n.lambda)
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
  
  W_inv.trn=diag(diag(S.trn)^(-0.5))
  
  R.trn=W_inv.trn%*%S.trn%*%W_inv.trn   #sample correlation matrix of training data
  
	 
	 for(j in 1:n.lambda)
	 {
	   lambda.mat=matrix(lambda[j],p,p)
     diag(lambda.mat)=0
     spice.trn=QUIC(R.trn,lambda.mat,msg = 0) #run spice
                   
       spice.icov.trn=W_inv.trn%*%spice.trn$X%*%W_inv.trn
  
          loss.lambda[j]=loss.lambda[j]+ norm( S.vld%*%spice.icov.trn - diag(1, p) ,"F")^2
          
          
	 }
 }


rm(X)
rm(X_vld)
rm(X_trn)
rm(S.trn)
rm(S.vld)
rm(spice.trn)
rm(spice.icov.trn)
rm(W_inv.trn)
rm(R.trn)


#optimal lambda
loss.lambda[which(loss.lambda==-Inf)]=NaN
index.lambda.opt=which(loss.lambda==min(loss.lambda,na.rm = TRUE)) 
lambda.opt=lambda[index.lambda.opt][1]

#spice for the whole data
lambda.opt.mat=matrix(lambda.opt,p,p)
diag(lambda.opt.mat)=0
spice=QUIC(R,lambda.opt.mat,msg = 0) #run spice

Omega.est=W_inv%*%spice$X%*%W_inv

delta=Omega-Omega.est

dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")

####### non-thresholding spice for Graphical model
##TPR/FPR
Omega.up=Omega[upper.tri(Omega,diag=F)]
Omega.est.up=Omega.est[upper.tri(Omega.est,diag=F)]


n.nonzero=sum(Omega.up!=0)
n.zero=sum(Omega.up==0)
TP=sum(Omega.est.up!=0 & Omega.up!=0)
FP=sum(Omega.est.up!=0 & Omega.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero

##output
time=proc.time()-ptm
output=c(dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR,n.nonzero,n.zero,time)

  

write.table(output,file=paste("model_band_50_F_0d01_spice_n",n,"_p",p,"_alpha",alpha,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)
