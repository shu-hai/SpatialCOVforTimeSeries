library("MASS")   #for mvrnorm
library("rwt") # for hard and soft thresholding
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) 


p=400#400

n=200



##############  Model 1 i.i.d. simulated data
Sigma=matrix(0,p,p) #precision matrix
Sigma=0.6^abs(row(Sigma)-col(Sigma)) #AR(1) matrix with parameter 0.6

## X
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) # seed for mvrnorm
X=mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma) #n*p matrix

X=t(X)

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


##################pertumation of X
X=X[,sample(1:n,n,replace=F)]

##################### 10-fold CV #######################

H_1=10


loss.hard=numeric(n.thresh)
loss.soft=numeric(n.thresh)



n.blk=n/H_1
X_i=(1:H_1)*n.blk



for(i in 1:H_1)
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

######################### 5-fold CV #######################

H_1=5


loss.hard=numeric(n.thresh)
loss.soft=numeric(n.thresh)



n.blk=n/H_1
X_i=(1:H_1)*n.blk



for(i in 1:H_1)
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



write.table(output,file=paste("model_ar_iid_per_n",n,"_p",p,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)



