######################################
### Load packages and source files ###
######################################
rm(list=ls())
library(parallel)
library(methods)
library(fields)
library(mvtnorm)
library(Rfast)
library(partitions)

wd <- "code/simulation/"
setwd(wd)

source('simu_Dombry_et_al.R',chdir=TRUE)
source('MLE_BrownResnick.R',chdir=TRUE)


###################################
### Main parameters (as inputs) ###
###################################

D <- 4^2 # dimension of the data
n <- 50 # sample size for numerical experiments
type <- 1 # type of correlation function (1=exponential, 2=powered exponential, 3=exchangeable)
par <- c(1,1) # parameter values
fixed <- c(FALSE,TRUE) # indicates whether or not the parameters are held fixed
estimator <- 2 # 1=MLE, 2=Composite likelihood estimator (MCLE), 3=Vecchia estimator (MVLE)
q.comp <- 2 # dimension of lower-dimensional margins used in the composite likelihoods
dmax <- 1.4143 # maximum distance for selecting composite likelihood components 
q.vecchia <- 1 # size of conditioning set in the Vecchia approximation (d=q-1)
ordering.vecchia <- 1 ## 1=vertical coordinate ordering, 2=random ordering, 3=middle-out ordering, 4=max-min ordering
R <- 10 # replications
ncores <- 10 # number of cores for parallel computing

set.seed(1)

x <- y <- seq(1,round(sqrt(D)),by=1)
loc <- as.matrix(expand.grid(x,y))
distmat <- rdist(loc)
rownames(distmat) <- colnames(distmat) <- NULL

##################
### Estimation ###
##################
# function for the r-th simulation experiment
mle.r <- function(r){
  
  ######################################
  ### Simulation of the r-th dataset ###
  ######################################
  vario <- function(x){
    sig <- Sigma(par,matrix(c(0,sqrt(sum(x^2)),sqrt(sum(x^2)),0),2,2),type)
    return(sig[1,1]-sig[1,2])
  }
  data <- simu_extrfcts(model="brownresnick",no.simu=n,coord=loc,vario=vario)$res # data simulation
  init <- par
  
  ####################################
	## Loop for Full likelihood (MLE) ##
	####################################
	if(estimator==1){
		tr <- try( time <- system.time( est <- MLE.BR(data=data,init=init,fixed=fixed,distmat=distmat,type=type)$par )[3] ) ### Full likelihood estimation
	  if(!is(tr,"try-error")){
	   est2 <- init
	   est2[!fixed] <- est
	   est <- est2
	  }
	}

  ##########################################
	## Loop for Composite likelihood (MCLE) ##
	##########################################
  if(estimator==2){
    all.index <- comb_n(D,q.comp) ## all possible q-uplets
    max.dist <- function(x){
      return(max(distmat[x,x]))
    }
    index <- all.index[,which(apply(all.index,2,max.dist)<=dmax)]
    
    tr <- try( time <- system.time( est <- MCLE.BR(data=data,init=init,fixed=fixed,distmat=distmat,type=type,index=index)$par )[3] ) ### Composite likelihood estimation based on q-uplets with maximum pairwise distance dmax
    if(!is(tr,"try-error")){
      est2 <- init
      est2[!fixed] <- est
      est <- est2
    }
  }
  
  ###########################################
  ## Loop for Vecchia approximation (MVLE) ##
  ###########################################
  if(estimator==3){
    vecchia.seq1<-c(1:D)[order(loc[,2])] ## Vecchia sequence based on vertical coordinate-ordering
    vecchia.seq2<-sample(1:D,D,replace=FALSE) ## Vecchia sequence based on completely random ordering
    vecchia.seq3<-order(distmat[which.min(apply(distmat,2,mean))[1],]) ## Vecchia sequence based on middle-out ordering
    vecchia.seq4<-c(which.min(distmat[which.min(colmeans(distmat))[1],])) ## Vecchia sequence based on max-min ordering
    for(j in 2:D){
      vecchia.seq4[j] <- c(1:D)[-vecchia.seq4][which.max(rowMins(matrix(distmat[-vecchia.seq4,vecchia.seq4],ncol=length(vecchia.seq4)),value=TRUE))[1]]
    }
    vecchia.seqs <- cbind(vecchia.seq1,vecchia.seq2,vecchia.seq3,vecchia.seq4)
    vecchia.seq <- vecchia.seqs[,ordering.vecchia]
    neighbours <- function(ind,vecchia.seq,q){
      if(ind==1){
        ind.neighbours <- rep(NA,q)
      }
      if(ind>=2){
        ind.ind <- vecchia.seq[ind] #index of ith-variable in the Vecchia sequence
        ind.past <- vecchia.seq[1:(ind-1)] #index of the "past" observations in the Vecchia sequence  
        d.past <- distmat[ind.past,vecchia.seq[ind]] #distance of the ith-variable to the "past" observations
        ind.neighbours <- ind.past[order(d.past)[1:q]] #choose "neighbours" as the closest observations in the "past"  
      }
      return(ind.neighbours)
    }
    neighbours.mat <- sapply(1:D,FUN=neighbours,vecchia.seq=vecchia.seq,q=q.vecchia)
    if(q.vecchia==1){
      neighbours.mat <- matrix(neighbours.mat,nrow=1)
    }
    
    tr <- try( time <- system.time( est <- MVLE.BR(data=data,init=init,fixed=fixed,distmat=distmat,type=type,vecchia.seq=vecchia.seq,neighbours=neighbours.mat)$par )[3] ) ### Vecchia likelihood estimation based on conditioning sets of size q.vecchia, using the lth-ordering
    if(!is(tr,"try-error")){
      est2 <- init
      est2[!fixed] <- est
      est <- est2
    }
  }
  
  return (list("estimates"=est,"time"=time))
}  

## Running the estimation (Warning: may take some time depending on the simulation setting!)
res <- mclapply(c(1:R),FUN=mle.r,mc.cores=ncores)

## Comparison between true parameter values and parameter estimates
errors <- matrix(nrow=R,ncol=length(par))
for(r in 1:R){
  errors[r,] <- res[[r]]$estimates-par
}
boxplot(errors)
abline(h=0,col="red")


