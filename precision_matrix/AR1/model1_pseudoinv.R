setwd("/Users/haishu/Desktop/work/project3_08072014/simulation/model1/result/pseudoinv")
library("MASS")   #for mvrnorm, ginv

#alpha=0.1
#p=100#400

alpha_r=c(0.1,0.25,0.5,1,2)
p_r=c(100,200,300,400)

n=200

for(alpha in alpha_r)
{
for(p in p_r)
{

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

for(seed in 1:100)
{
#seed=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
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

##sample covariance matrix
S=cov(t(X))*(n-1)/n

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


write.table(output,file=paste("model1_pseudoinv_n",n,"_p",p,"_alpha",alpha,"_seed",seed,".txt",sep=""),row.names = FALSE,col.names = FALSE)
}
}
}
