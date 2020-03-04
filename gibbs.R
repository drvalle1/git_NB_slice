rm(list=ls(all=TRUE)) 
library('Rcpp')
set.seed(1)

#get functions
setwd('U:\\GIT_models\\git_NB_slice')
sourceCpp('slice_betas.cpp')
sourceCpp('slice_NBN.cpp')

#get data
dat=read.csv('fake data.csv',as.is=T)
ind=grep('y',colnames(dat))
y=data.matrix(dat[,ind])
xmat=data.matrix(cbind(1,dat[,-ind]))

#get initial values
ngroup=4
nparam=ncol(xmat)
betas=matrix(0,nparam,ngroup)
NBN=10
nloc=nrow(y)

#basic settings
ngibbs=1000
nburn=ngibbs/2

#priors
var.betas=c(10,rep(1,ncol(xmat)-1))

#to store outcomes from gibbs sampler
betas.out=matrix(NA,ngibbs,nparam*ngroup)
NBN.out=matrix(NA,ngibbs,1)

#useful stuff for slice sampler algorithm
w.betas=0.1
w.NBN=10

#run gibbs sampler
options(warn=2)
for (i in 1:ngibbs){
  print(i)
  
  #sample betas
  betas=SampleBetas(param=betas,y=y,xmat=xmat,w=w.betas,nparam=nparam,ncomm=ngroup,var1=var.betas,NBN=NBN)
  
  #sample NBN
  media=exp(xmat%*%betas) #get mean
  NBN=SampleNBN(Media=media,y=y,NBN=NBN,w=w.NBN)
  
  #store results  
  betas.out[i,]=betas
  NBN.out[i]=NBN
}

plot(NBN.out,type='l')
quantile(NBN.out,c(0.025,0.5,0.975))

rango=range(c(betas.true,betas))
plot(betas.true,betas,xlim=rango,ylim=rango)
lines(rango,rango,col='red')