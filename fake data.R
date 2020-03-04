rm(list=ls(all=TRUE)) 
set.seed(2)

n=100 #bigger n = smaller variance
nobs=1000
ncomm=4
nparam=4
xmat=matrix(runif(ncomm*((nparam)-1),min=-1,max=1),nobs,nparam-1)
colnames(xmat)=paste0('cov',1:(nparam-1))
xmat1=cbind(1,xmat)
betas.true=param=matrix(runif(nparam*ncomm,max=2),nparam,ncomm)
mean1=exp(xmat1%*%param); range(mean1)

tmp=rnbinom(nobs*ncomm,mu=mean1,size=n)
y=matrix(tmp,nobs,ncomm)
colnames(y)=paste0('y',1:ncomm)

fim=cbind(y,xmat)

setwd('U:\\GIT_models\\git_NB_slice')
write.csv(fim,'fake data.csv',row.names=F)