args <- commandArgs(TRUE)
load("data/data.RData");source("code/MLE_BrownResnick.R"); 
#default setting
q.comp <- 2 # number of nearest neigbours for composite likelihood c(2,3)
dmax = 2
q.vecchia <- 2 # number of 'historical sites' in the vecchia method 1,2,3,4
vecchia.order = 1 # vecchia ordering c(1,2,3,4)
ratio=1 #ratio of the dataset locations used for training, 1,1/2,1/4,1/8,1/16.
init = c(1,100,1,0)
fixed = c(F,F,T,T)
for (arg in args) eval(parse(text = arg))

### Libraries that needed for the functions ###
library(evd)
library(parallel)
library(methods)
library(mvtnorm)
library(partitions)
ncores = detectCores()

## select the locations according to max-min ordering accoording to the ratio 
distmat <- as.matrix(dist(coord))
loc.sub.trans = coord
ind.sub <- c()
ind.sub[1]<-c(which.min(distmat[which.min(colMeans(distmat))[1],])) ## Vecchia sequence based on max-min ordering
for(j in 2:nrow(loc.sub.trans)){
  ind.sub[j] <- c(1:nrow(loc.sub.trans))[-ind.sub[1:(j-1)] ][which.max(
		apply(matrix(distmat[-ind.sub[1:(j-1)],ind.sub[1:(j-1)]],ncol=(j-1)),1,min))[1]]
}
if(ratio>0.5){
  ind.sub <- rev(ind.sub)[1:floor(nrow(loc.sub.trans)*ratio)]
}else{
  ind.sub <- ind.sub[1:floor(nrow(loc.sub.trans)*ratio)]
}


set.seed(324234)
## fit the model using composite likelihood
D <- length(ind.sub)
distmat.temp <- distmat[ind.sub,ind.sub]
all.index <- combn(x=D,m=q.comp)
max.dist <- function(id){
    x = all.index[,id]
    return(max(distmat.temp[x,x]))
}
dist.max <- order(unlist(mclapply(1:ncol(all.index),max.dist,mc.cores = ncores)))  
dist.ind<-which(dist.max <= D*dmax)
all.index = all.index[,dist.ind]
time.used <- system.time( fit.result <- MCLE.BR(data=samples.skew.normal1[[1]][,ind.sub],init=init,fixed=fixed,distmat=loc.sub.trans[ind.sub,],FUN = vario.func,index=all.index,ncores,method="Brent",maxit=1000,hessian=FALSE))[3]


## fit the model using Vecchia method
### Vecchia approximation with 4 different Vecchia sequences and 3 different Vecchia orders (Brown-Resnick with exponential correlation)
D <- length(ind.sub)
distmat.temp <- distmat[ind.sub,ind.sub]
vecchia.seq <- matrix(NA,ncol=4,nrow=D)
vecchia.seq[,1]<-c(1:D)[order(loc.sub.trans[ind.sub,2])] ## Vecchia sequence based on vertical coordinate-ordering
vecchia.seq[,2]<-sample(1:D,D,replace=FALSE) ## Vecchia sequence based on completely random ordering
vecchia.seq[,3]<-order(distmat.temp[which.min(apply(distmat.temp,2,mean))[1],]) ## Vecchia sequence based on middle-out ordering
vecchia.seq[1,4]<-c(which.min(distmat.temp[which.min(colMeans(distmat.temp))[1],])) ## Vecchia sequence based on max-min ordering
for(j in 2:D){
    vecchia.seq[j,4] <- c(1:D)[-vecchia.seq[1:(j-1),4] ][which.max(apply(matrix(distmat.temp[-vecchia.seq[1:(j-1),4],vecchia.seq[1:(j-1),4]],ncol=(j-1)),1,min))[1]]
}

neighbours.mat <- sapply(1:D,FUN=neighbours,vecchia.seq=vecchia.seq[,vecchia.order],
					q=q.vecchia,distmat=distmat.temp)
if(q.vecchia==1){
    neighbours.mat <- matrix(neighbours.mat,nrow=1)
}
time.used <- system.time(fit.result <- MVLE.BR(data=maxima.frechet[,ind.sub],init=init,fixed=fixed,distmat=loc.sub.trans[ind.sub,],FUN = vario.func,vecchia.seq=vecchia.seq[,vecchia.order],neighbours=neighbours.mat,ncores,method="Nelder-Mead",maxit=1000,hessian=hessian) )[3]

## calculate the prediction power
nlog.pred.val = NULL
if(model.eval & ratio!=1){
  set.seed(1)
  neighbours <- function(ind,ind.sub,qr){
    ind.neighbours <- ind.sub[order(distmat[ind,ind.sub])[1:qr]]
    return(ind.neighbours)
  }
  neighbours2predict <- sapply((1:nrow(distmat))[-ind.sub],FUN=neighbours,ind.sub=ind.sub,qr=4)
  par.val <- fit.result$par
  nlog.pred.val <- nlogVecchialik.BR(par=par.val,data=maxima.frechet,distmat=loc.sub.trans,FUN=vario.func,vecchia.seq=(1:nrow(distmat))[-ind.sub],neighbours = neighbours2predict,ncores=ncores)
}





