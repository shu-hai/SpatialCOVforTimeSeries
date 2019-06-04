setwd("/Users/haishu/Desktop/work/project3_08072014/simulation/model1/result/pseudoinv")
library("MASS")   #for mvrnorm, ginv

p_r=c(100,200,300,400)
for(p in p_r)
{
#p=100#400

n=200

##############  Model 1 i.i.d. simulated data
Omega=matrix(0,p,p) #precision matrix
Omega=0.6^abs(row(Omega)-col(Omega)) #AR(1) matrix with parameter 0.6

for(seed in 1:100)
{
## X
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


output=c(dist.norm.1,dist.norm.2,dist.norm.F)


write.table(output,file=paste("model1_iid_pseudoinv_n",n,"_p",p,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)
}
}
