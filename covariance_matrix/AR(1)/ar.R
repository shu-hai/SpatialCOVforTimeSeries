library("MASS")   #for mvrnorm
library("rwt") # for hard and soft thresholding
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) 


alpha=0.1
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
Sigma=matrix(0,p,p) #precision matrix
Sigma=0.6^abs(row(Sigma)-col(Sigma)) #AR(1) matrix with parameter 0.6

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
corr=matrix(0,p,p)
corr=0.6^abs(row(corr)-col(corr)) #AR(1) matrix with parameter 0.6

#sample correlation
corr.hat=cor(t(X))

delta=corr-corr.hat

dist.norm.2=base::norm(delta,"2")
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
output=c(dist.norm.2,dist.norm.F,dist.norm.1)


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
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

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
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

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
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

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
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

delta=corr-corr.soft
dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
output=c(output,dist.norm.2,dist.norm.F,dist.norm.1)

write.table(output,file=paste("model_ar_n",n,"_p",p,"_alpha",alpha,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)



