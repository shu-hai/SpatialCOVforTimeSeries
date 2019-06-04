library("MASS")   #for mvrnorm, ginv
setwd("/Users/haishu/Desktop/work/project3_08072014/simulation/band_model/result/pseudoinv")

p_r=c(100,200,300,400)
for(p in p_r)
{
n=200

##############  Model band i.i.d. simulated data
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


## X
for(seed in 1:100)
{
#seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(seed) # seed for mvrnorm
X=mvrnorm(n=n,mu=rep(0,p),Sigma=solve(Omega)) #n*p matrix

##sample covariance matrix
S=cov(X)*(n-1)/n

## pseudo-inverse S
Omega.est=tryCatch({
  ginv(S)  
}, error = function(e){
  tryCatch({
    ginv(S/sd(S))/sd(S)#scale the matrix to make ginv stable
  }, error = function(e){  
    tryCatch({
      ginv(S/max(abs(S)))/max(abs(S))
    },error = function(e){       
      
      ginv(S/mean(diag(S)))/mean(diag(S))
      
    })
  })
})

delta=Omega-Omega.est

dist.norm.1=norm(delta,"1")
dist.norm.F=norm(delta,"F")
dist.norm.2=base::norm(delta,"2")
##TPR/FPR
Omega.up=Omega[upper.tri(Omega,diag=F)]
Omega.est.up=Omega.est[upper.tri(Omega.est,diag=F)]


n.nonzero=sum(Omega.up!=0)
n.zero=sum(Omega.up==0)
TP=sum(Omega.est.up!=0 & Omega.up!=0)
FP=sum(Omega.est.up!=0 & Omega.up==0)

TPR=TP/n.nonzero
FPR=FP/n.zero

output=c(dist.norm.1,dist.norm.2,dist.norm.F,TPR,FPR,n.nonzero,n.zero)

write.table(output,file=paste("model_band_iid_pseudoinv_n",n,"_p",p,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)

}
}
