library("MASS")   #for mvrnorm
library("rwt") # for hard and soft thresholding
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) 


p=300#400

n=200

##############  Model band i.i.d. simulated data
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



## X
seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) # seed for mvrnorm
X=mvrnorm(n=n,mu=rep(0,p),Sigma=Sigma) #n*p matrix

X=t(X)

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
thresh=seq(0.01,0.99,by=0.01)
n.thresh=length(thresh)

##################pertumation of X
X=X[,sample(1:n,n,replace=F)]


######################### 10-fold CV  #######################

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


write.table(output,file=paste("model_mr_iid_per_n",n,"_p",p,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)



