library("MASS")   #for mvrnorm
library("rwt") # for hard and soft thresholding
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) 



alpha=2
p=300#400

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
####construct Sigma matrix
Sigma=matrix(0,p,p)
for(i in 1:(p-2))
{
    
    Sigma[i,i]=0.5
    Sigma[i,(i+1)]=0.6
    Sigma[i,(i+2)]=0.3
    
}
Sigma[p-1,p-1]=0.5
Sigma[p,p]=0.5
Sigma[p-1,p]=0.6
Sigma=Sigma+t(Sigma)




####generate Y and e vector
num=N*n
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) # seed for mvrnorm
mvn.vec=mvrnorm(n=num,mu=rep(0,p),Sigma=Sigma) #num*p matrix

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



#the true correlation matrix
corr=Sigma


#sample correlation
corr.hat=cor(t(X))

delta=corr-corr.hat

dist.norm.2=base::norm(delta,"2")
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")

#
corr.up=corr[upper.tri(corr,diag=F)]
n.nonzero=sum(corr.up!=0)
n.zero=sum(corr.up==0)


corr.hat.up=corr.hat[upper.tri(corr.hat,diag=F)]

TP=sum(corr.hat.up!=0 & corr.up!=0)
FP=sum(corr.hat.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
#
output=c(dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)


##threshold parameter
thresh=seq(0.01,0.5,by=0.01)
n.thresh=length(thresh)

############ H_1=H_2=10, intuitive CV ###############
H_1=10
H_2=10

loss.hard=numeric(n.thresh)
loss.soft=numeric(n.thresh)



n.blk=n/H_1
X_i=c((1:H_1)*n.blk,sample(n.blk:n,H_2))

H.12=H_1+H_2

for(i in 1:H.12)
{
  X_vld_e=X_i[i]
  X_vld_b=X_i[i]-(n.blk-1)
  
  X_vld=X[,X_vld_b:X_vld_e]# validation data
  
  X_trn_l=X_vld_b-(n.blk+1)
  X_trn_r=X_vld_e+(n.blk+1)
  
  if(X_trn_l<1){
    X_trn=X[,X_trn_r:n] 
  }else if(X_trn_r>n){
    X_trn=X[,1:X_trn_l]
  }else{
    X_trn=cbind(X[,1:X_trn_l],X[,X_trn_r:n])# training data
  }
  
  corr.vld=cor(t(X_vld)) 
  corr.trn=cor(t(X_trn))
  
  for(j in 1:n.thresh)
  {
    corr.trn.hard=hardTh(corr.trn,thresh[j])
    corr.trn.soft=softTh(corr.trn,thresh[j])
    diag(corr.trn.soft)=1
    
    loss.hard[j]=loss.hard[j]+norm(corr.trn.hard-corr.vld,"F")^2
    loss.soft[j]=loss.soft[j]+norm(corr.trn.soft-corr.vld,"F")^2
  }
  
}

index.hard.opt=which(loss.hard==min(loss.hard))
index.soft.opt=which(loss.soft==min(loss.soft))

t.hard.opt=thresh[index.hard.opt[1]]
t.soft.opt=thresh[index.soft.opt[1]]



corr.hard=hardTh(corr.hat,t.hard.opt)
corr.soft=softTh(corr.hat,t.soft.opt)
diag(corr.soft)=1


delta=corr-corr.hard
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")

corr.hard.up=corr.hard[upper.tri(corr.hard,diag=F)]

TP=sum(corr.hard.up!=0 & corr.up!=0)
FP=sum(corr.hard.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)

delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")

corr.soft.up=corr.soft[upper.tri(corr.soft,diag=F)]

TP=sum(corr.soft.up!=0 & corr.up!=0)
FP=sum(corr.soft.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)

############ H_1=5,H_2=10, intuitive CV ###############
H_1=5
H_2=10

loss.hard=numeric(n.thresh)
loss.soft=numeric(n.thresh)



n.blk=n/H_1
X_i=c((1:H_1)*n.blk,sample(n.blk:n,H_2))

H.12=H_1+H_2

for(i in 1:H.12)
{
  X_vld_e=X_i[i]
  X_vld_b=X_i[i]-(n.blk-1)
  
  X_vld=X[,X_vld_b:X_vld_e]# validation data
  
  X_trn_l=X_vld_b-(n.blk+1)
  X_trn_r=X_vld_e+(n.blk+1)
  
  if(X_trn_l<1){
    X_trn=X[,X_trn_r:n] 
  }else if(X_trn_r>n){
    X_trn=X[,1:X_trn_l]
  }else{
    X_trn=cbind(X[,1:X_trn_l],X[,X_trn_r:n])# training data
  }
  
  corr.vld=cor(t(X_vld)) 
  corr.trn=cor(t(X_trn))
  
  for(j in 1:n.thresh)
  {
    corr.trn.hard=hardTh(corr.trn,thresh[j])
    corr.trn.soft=softTh(corr.trn,thresh[j])
    diag(corr.trn.soft)=1
    
    loss.hard[j]=loss.hard[j]+norm(corr.trn.hard-corr.vld,"F")^2
    loss.soft[j]=loss.soft[j]+norm(corr.trn.soft-corr.vld,"F")^2
  }
  
}

index.hard.opt=which(loss.hard==min(loss.hard))
index.soft.opt=which(loss.soft==min(loss.soft))

t.hard.opt=thresh[index.hard.opt[1]]
t.soft.opt=thresh[index.soft.opt[1]]



corr.hard=hardTh(corr.hat,t.hard.opt)
corr.soft=softTh(corr.hat,t.soft.opt)
diag(corr.soft)=1


delta=corr-corr.hard
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
corr.hard.up=corr.hard[upper.tri(corr.hard,diag=F)]

TP=sum(corr.hard.up!=0 & corr.up!=0)
FP=sum(corr.hard.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)



delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
corr.soft.up=corr.soft[upper.tri(corr.soft,diag=F)]

TP=sum(corr.soft.up!=0 & corr.up!=0)
FP=sum(corr.soft.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)

######################### non-block H_1=H_2=10 ordinary CV #######################

H_1=10
H_2=10

loss.hard=numeric(n.thresh)
loss.soft=numeric(n.thresh)



n.blk=n/H_1
X_i=c((1:H_1)*n.blk,sample(n.blk:n,H_2))

H.12=H_1+H_2

for(i in 1:H.12)
{
  X_vld_e=X_i[i]
  X_vld_b=X_i[i]-(n.blk-1)
  
  X_vld=X[,X_vld_b:X_vld_e]# validation data
  
  X_trn_l=X_vld_b-1
  X_trn_r=X_vld_e+1
  
  if(X_trn_l<1){
    X_trn=X[,X_trn_r:n] 
  }else if(X_trn_r>n){
    X_trn=X[,1:X_trn_l]
  }else{
    X_trn=cbind(X[,1:X_trn_l],X[,X_trn_r:n])# training data
  }
  
  corr.vld=cor(t(X_vld)) 
  corr.trn=cor(t(X_trn))
  
  for(j in 1:n.thresh)
  {
    corr.trn.hard=hardTh(corr.trn,thresh[j])
    corr.trn.soft=softTh(corr.trn,thresh[j])
    diag(corr.trn.soft)=1
    
    loss.hard[j]=loss.hard[j]+norm(corr.trn.hard-corr.vld,"F")^2
    loss.soft[j]=loss.soft[j]+norm(corr.trn.soft-corr.vld,"F")^2
  }
  
}

index.hard.opt=which(loss.hard==min(loss.hard))
index.soft.opt=which(loss.soft==min(loss.soft))

t.hard.opt=thresh[index.hard.opt[1]]
t.soft.opt=thresh[index.soft.opt[1]]



corr.hard=hardTh(corr.hat,t.hard.opt)
corr.soft=softTh(corr.hat,t.soft.opt)
diag(corr.soft)=1


delta=corr-corr.hard
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
corr.hard.up=corr.hard[upper.tri(corr.hard,diag=F)]

TP=sum(corr.hard.up!=0 & corr.up!=0)
FP=sum(corr.hard.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)



delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
corr.soft.up=corr.soft[upper.tri(corr.soft,diag=F)]

TP=sum(corr.soft.up!=0 & corr.up!=0)
FP=sum(corr.soft.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)

######################### non-block H_1=5, H_2=10 ordinary CV #######################

H_1=5
H_2=10

loss.hard=numeric(n.thresh)
loss.soft=numeric(n.thresh)



n.blk=n/H_1
X_i=c((1:H_1)*n.blk,sample(n.blk:n,H_2))

H.12=H_1+H_2

for(i in 1:H.12)
{
  X_vld_e=X_i[i]
  X_vld_b=X_i[i]-(n.blk-1)
  
  X_vld=X[,X_vld_b:X_vld_e]# validation data
  
  X_trn_l=X_vld_b-1
  X_trn_r=X_vld_e+1
  
  if(X_trn_l<1){
    X_trn=X[,X_trn_r:n] 
  }else if(X_trn_r>n){
    X_trn=X[,1:X_trn_l]
  }else{
    X_trn=cbind(X[,1:X_trn_l],X[,X_trn_r:n])# training data
  }
  
  corr.vld=cor(t(X_vld)) 
  corr.trn=cor(t(X_trn))
  
  for(j in 1:n.thresh)
  {
    corr.trn.hard=hardTh(corr.trn,thresh[j])
    corr.trn.soft=softTh(corr.trn,thresh[j])
    diag(corr.trn.soft)=1
    
    loss.hard[j]=loss.hard[j]+norm(corr.trn.hard-corr.vld,"F")^2
    loss.soft[j]=loss.soft[j]+norm(corr.trn.soft-corr.vld,"F")^2
  }
  
}

index.hard.opt=which(loss.hard==min(loss.hard))
index.soft.opt=which(loss.soft==min(loss.soft))

t.hard.opt=thresh[index.hard.opt[1]]
t.soft.opt=thresh[index.soft.opt[1]]



corr.hard=hardTh(corr.hat,t.hard.opt)
corr.soft=softTh(corr.hat,t.soft.opt)
diag(corr.soft)=1


delta=corr-corr.hard
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
corr.hard.up=corr.hard[upper.tri(corr.hard,diag=F)]

TP=sum(corr.hard.up!=0 & corr.up!=0)
FP=sum(corr.hard.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR)



delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
corr.soft.up=corr.soft[upper.tri(corr.soft,diag=F)]

TP=sum(corr.soft.up!=0 & corr.up!=0)
FP=sum(corr.soft.up!=0 & corr.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1,TPR,FPR,n.nonzero,n.zero)

write.table(output,file=paste("model_mr_n",n,"_p",p,"_alpha",alpha,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)



