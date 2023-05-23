#! /usr/bin/Rscript
args <- commandArgs(TRUE)
for (arg in args) eval(parse(text = arg))
rm(arg, args)

######################################
### Load packages and source files ###
######################################

library(parallel)
library(methods)
library(fields)
library(mvtnorm)
library(Rfast)
library(partitions)

wd <- "."
setwd(wd)

source('code/simu_Dombry_et_al.R')
source('code/MLE_BrownResnick.R')


###################################
### Main parameters (as inputs) ###
###################################

Ds <- 10^2 # dimension of the data
ns <- c(100) # sample size for numerical experiments
alpha = 1;lambda = c(1,2);a=0.5;theta= pi/8
pars <- t(as.matrix(expand.grid(alpha,lambda,a,theta))) # Different parameter values
fixed <- c(FALSE,FALSE,FALSE,FALSE) # indicates whether or not the parameters are held fixed
estimators <- 3 # 1=MLE, 2=Composite likelihood estimator (MCLE), 3=Vecchia estimator (MVLE)
qs.comp <- c(2,3,4,5) # dimension of lower-dimensional margins used in the composite likelihoods
dmaxs <- c(1.0001,1.4143,2.0001) ## maximum distances for selecting composite likelihood components (here, from 1st to 5th-order neighbors when Dx=Dy=1)
qs.vecchia <- c(2,5) #sizes of conditioning set in the Vecchia approximation
orderings.vecchia <- c(1,2,4) ## 1=vertical coordinate ordering, 2=random ordering, 3=middle-out ordering, 4=max-min ordering
MDA <- FALSE # should the data be simulated in the max-domain of attraction?
M <- 100 # block size, if simulation is done in the MDA
Rs <- c(1:20)*20 # replications
ncores <- 20 # number of cores for parallel computing
simul <- 1 
subsimul <- 'a'

### Writing down the inputs into a file to remember it
if(any(Rs==1)){ ## This prevents conflicts with parallel computing...
  parameters <- list(Ds,ns,pars,fixed,qs.comp,dmaxs,orderings.vecchia,qs.vecchia)
  names(parameters) <- c("Ds","ns","pars","fixed","qs.comp","dmaxs","orderings.vecchia","qs.vecchia")
   for(i in 1:length(parameters)){
    out <- paste(names(parameters)[[i]]," <- c(",paste(parameters[[i]],collapse= ","),")\n",sep="")
    cat(out,file=paste("data/Simulation",simul,"_SimulationParameters.txt",sep=""),append=TRUE)
  } 
}

### Putting parameter settings in matrix form, and extracting number of parameter settings and dimensionality of parameter vector
if(is.matrix(pars)){
  npars <- ncol(pars) ## number of parameter settings
  dimpars <- nrow(pars) ## dimensionality of parameter vector
} else{
  npars <- length(pars) ## number of parameter settings
  dimpars <- 1 ## dimensionality of parameter vector
}
pars.mat <- matrix(pars,nrow=dimpars,ncol=npars)

###########################################
### Simulation study for Gaussian model ###
###########################################

#Simulation for a single replicate r in input vector Rs
mle.r <- function(r,ncores=1){
  
  if(any(estimators==1)){
    estimates.MLE.r <- array(dim=c(length(Ds),length(ns),npars,dimpars))
	  times.MLE.r <- array(dim=c(length(Ds),length(ns),npars))
  }
  if(any(estimators==2)){
    estimates.Comp.r <- array(dim=c(length(Ds),length(ns),npars,length(qs.comp),length(dmaxs),dimpars))
	  times.Comp.r <- array(dim=c(length(Ds),length(ns),npars,length(qs.comp),length(dmaxs)))
  }
  if(any(estimators==3)){
    estimates.Vecchia.r <- array(dim=c(length(Ds),length(ns),npars,length(qs.vecchia),length(orderings.vecchia),dimpars))
	  times.Vecchia.r <- array(dim=c(length(Ds),length(ns),npars,length(qs.vecchia),length(orderings.vecchia)))
  }
  
	for(i in 1:length(Ds)){
		D <- Ds[i]
		x <- y <- seq(1,round(sqrt(D)),by=1)
		loc <- as.matrix(expand.grid(x,y))
		distmat <- rdist(loc)
		#distmat <- as.matrix(dist(loc))
		rownames(distmat) <- colnames(distmat) <- NULL
		D <- nrow(loc)
		
		####################################
		## Loop for Full likelihood (MLE) ##
		####################################
		if(any(estimators==1)){
		  for(j in 1:length(ns)){
		    n <- ns[j]
		   	for(k in 1:npars){
		      par <- pars.mat[,k]
		      set.seed(18462*r+8934+r) ## we fix the random seed for each D, par, n, but we take a different random seed for each experiment r=1,...,1024 (so the Vecchia sequence based on random ordering changes for each r)
		      if(!MDA){
		        vario <- function(x){ ## coordinates
		          if(!is.matrix(x)){
		            sig <- vario.func(loc=matrix(x,ncol=2),init,relocate = FALSE)
		            val=sig[1,1]/2
		          }else{
		            sig <- vario.func(loc=x,init,relocate = FALSE)  
		            val = sig[1,1]/2+sig[2,2]/2-sig[1,2]
		          }
		          return(val)
		        }
		        data <- simu_extrfcts(model="brownresnick",no.simu=n,coord=loc,vario=vario)$res # data simulation
		      } else{
		        N <- M*n
		        sigma2 <- drop(vario.func(loc=c(0,0),par,relocate = FALSE))
		        lambda <- sqrt(sigma2)
		        W <- mvtnorm::rmvnorm(N,sigma=cov2cor(vario.func(loc,par,relocate = FALSE))) + rexp(N,rate=lambda) # data simulation  from factor copula model
		        Copula <- pnorm(W)-exp(lambda^2/2-lambda*W)*pnorm(W-lambda)
		        data.MDA <- 1/(1-Copula)
		        data <- matrix(nrow=n,ncol=D) # this will contain block maxima
		        for(b in 1:n){
		          data[b,] <- colMaxs(data.MDA[((b-1)*M+1):(b*M),],value=TRUE)/M #componentwise maxima for block b.
		          #data[b,] <- apply(data.MDA[((b-1)*M+1):(b*M),],2,max)/M #componentwise maxima for block b.
		        }
		      }
		      
		      init <- par; init[!is.finite(init)] = 0
		      
		      tr <- try( time <- system.time( est <- MLE.BR(data=data,init=init,fixed=fixed,distmat=loc,FUN=vario.func)$par )[3] ) ### Full likelihood estimation
		      if(!is(tr,"try-error")){
		        est2 <- init
		        est2[!fixed] <- est
		        estimates.MLE.r[i,j,k,] <- est2
		        times.MLE.r[i,j,k] <- time
		      }
		    }
		  }
		}
		##########################################
		## Loop for Composite likelihood (MCLE) ##
		##########################################
		if(any(estimators==2)){
		  for(l in 1:length(qs.comp)){
		    q.comp <- qs.comp[l]
		    all.index <- comb_n(D,q.comp) ## all possible q-uplets
		    max.dist <- function(x){
		      return(max(distmat[x,x]))
		    }
		    max.dist <- apply(all.index,2,max.dist) ## maximum pairwise distance for all q-uplets
		    for(m in 1:length(dmaxs)){
		      dmax <- dmaxs[m]
		      index <- all.index[,which(max.dist<=dmax)]
		      for(j in 1:length(ns)){
		        n <- ns[j]
		        for(k in 1:npars){
		          par <- pars.mat[,k]
		          set.seed(18462*r+8934+r) ## we fix the random seed for each D, par, n, but we take a different random seed for each experiment r=1,...,1024 (so the Vecchia sequence based on random ordering changes for each r)
		          if(!MDA){
		            vario <- function(x){ ## coordinates
		              if(!is.matrix(x)){
		                sig <- vario.func(loc=matrix(x,ncol=2),init,relocate = FALSE)
		                val=sig[1,1]/2
		              }else{
		                sig <- vario.func(loc=x,init,relocate = FALSE)  
		                val = sig[1,1]/2+sig[2,2]/2-sig[1,2]
		              }
		              return(val)
		            }
		            data <- simu_extrfcts(model="brownresnick",no.simu=n,coord=loc,vario=vario)$res # data simulation
		          } else{
		            N <- M*n
		            sigma2 <- drop(vario.func(loc=c(0,0),par,relocate = FALSE))
		            lambda <- sqrt(sigma2)
		            W <- mvtnorm::rmvnorm(N,sigma=cov2cor(vario.func(loc,par,relocate = FALSE))) + rexp(N,rate=lambda) # data simulation from factor copula model
		            Copula <- pnorm(W)-exp(lambda^2/2-lambda*W)*pnorm(W-lambda)
		            data.MDA <- 1/(1-Copula)
		            data <- matrix(nrow=n,ncol=D) # this will contain block maxima
		            for(b in 1:n){
		              data[b,] <- colMaxs(data.MDA[((b-1)*M+1):(b*M),],value=TRUE)/M #componentwise maxima for block b.
		              #data[b,] <- apply(data.MDA[((b-1)*M+1):(b*M),],2,max)/M #componentwise maxima for block b.
		            }
		          }
		          
		          init <- par; init[!is.finite(init)] = 0
		      
		          tr <- try( time <- system.time( est <- MCLE.BR(data=data,init=init,fixed=fixed,distmat=loc,FUN=vario.func,index=index,ncores=ncores)$par )[3] ) ### Composite likelihood estimation based on q-uplets with maximum pairwise distance dmax
		          if(!is(tr,"try-error")){
		            est2 <- init
		            est2[!fixed] <- est
		            estimates.Comp.r[i,j,k,l,m,] <- est2
		            times.Comp.r[i,j,k,l,m] <- time
		          }
		        }
		      }
		    }
		    rm(all.index)
		  }
		}
		
		###########################################
		## Loop for Vecchia approximation (MVLE) ##
		###########################################
		if(any(estimators==3)){
		  vecchia.seq1<-c(1:D)[order(loc[,2])] ## Vecchia sequence based on vertical coordinate-ordering
		  vecchia.seq2<-sample(1:D,D,replace=FALSE) ## Vecchia sequence based on completely random ordering
		  vecchia.seq3<-order(distmat[which.min(apply(distmat,2,mean))[1],]) ## Vecchia sequence based on middle-out ordering
		  vecchia.seq4<-c(which.min(distmat[which.min(colmeans(distmat))[1],])) ## Vecchia sequence based on max-min ordering
		  for(j in 2:D){
		    vecchia.seq4[j] <- c(1:D)[-vecchia.seq4][which.max(rowMins(matrix(distmat[-vecchia.seq4,vecchia.seq4],ncol=length(vecchia.seq4)),value=TRUE))[1]]
		  }
		  vecchia.seqs <- cbind(vecchia.seq1,vecchia.seq2,vecchia.seq3,vecchia.seq4)
		  for(l in 1:length(qs.vecchia)){
		    q.vecchia <- qs.vecchia[l]
		    for(m in 1:length(orderings.vecchia)){
		      vecchia.seq <- vecchia.seqs[,orderings.vecchia[m]]
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
		      for(j in 1:length(ns)){
		        n <- ns[j]
		        for(k in 1:npars){
		          par <- pars.mat[,k]
		          set.seed(18462*r+8934+r) ## we fix the random seed for each D, alpha, n, but we take a different random seed for each experiment r=1,...,1024 (so the Vecchia sequence based on random ordering changes for each r)
		          if(!MDA){
		            vario <- function(x){
		              sig <- vario.func(loc=rbind(x,0),par,relocate=FALSE)
		              return(sig[1,1]-sig[1,2])
		            }
		            data <- simu_extrfcts(model="brownresnick",no.simu=n,coord=loc,vario=vario)$res # data simulation
		          } else{
		            N <- M*n
		            sigma2 <- drop(vario.func(loc = c(0,0),par,relocate = FALSE))
		            lambda <- sqrt(sigma2)
		            W <- mvtnorm::rmvnorm(N,sigma=cov2cor(vario.func(loc,par,relocate = FALSE))) + rexp(N,rate=lambda) # data simulation from factor copula model
		            Copula <- pnorm(W)-exp(lambda^2/2-lambda*W)*pnorm(W-lambda)
		            data.MDA <- 1/(1-Copula)
		            data <- matrix(nrow=n,ncol=D) # this will contain block maxima
		            for(b in 1:n){
		              data[b,] <- colMaxs(data.MDA[((b-1)*M+1):(b*M),],value=TRUE)/M #componentwise maxima for block b.
		              #data[b,] <- apply(data.MDA[((b-1)*M+1):(b*M),],2,max)/M #componentwise maxima for block b.
		            }
		          }
		          
		          init <- par; init[!is.finite(init)] = 0
		      
		          tr <- try( time <- system.time( est <- MVLE.BR(data=data,init=init,fixed=fixed,distmat=loc,FUN = vario.func,vecchia.seq=vecchia.seq,neighbours=neighbours.mat,ncores = ncores)$par )[3] ) ### Vecchia likelihood estimation based on conditioning sets of size q.vecchia, using the lth-ordering
		          if(!is(tr,"try-error")){
		            est2 <- init
		            est2[!fixed] <- est
		            estimates.Vecchia.r[i,j,k,l,m,] <- est2
		            times.Vecchia.r[i,j,k,l,m] <- time
		          }
		        }
		      }
		    }
		  }
		}
	}	
	
  res <- list()
  if(any(estimators==1)){
	  res$estimates.MLE.r <- estimates.MLE.r
	  res$times.MLE.r <- times.MLE.r
  }
  if(any(estimators==2)){
    res$estimates.Comp.r <- estimates.Comp.r
    res$times.Comp.r <- times.Comp.r
  }
  if(any(estimators==3)){
    res$estimates.Vecchia.r <- estimates.Vecchia.r
    res$times.Vecchia.r <- times.Vecchia.r
  }
	  
	return(res)
}


results.list <- mclapply(Rs,FUN=mle.r,mc.cores=ncores)

 
#################################
### Organize and save results ###
#################################

if(any(estimators==1)){
  estimates.MLE <- array(dim=c(length(Rs),length(Ds),length(ns),npars,dimpars))
  times.MLE <- array(dim=c(length(Rs),length(Ds),length(ns),npars))
  dimnames(estimates.MLE)<-list(paste('R=',Rs,sep=''),paste('D=',Ds,sep=''),paste('n=',ns,sep=''),paste('par=c(',apply(pars.mat,2,paste,collapse=','),')',sep=''),paste('par ',c(1:dimpars),sep=''))
  dimnames(times.MLE)<-list(paste('R=',Rs,sep=''),paste('D=',Ds,sep=''),paste('n=',ns,sep=''),paste('par=c(',apply(pars.mat,2,paste,collapse=','),')',sep=''))
}
if(any(estimators==2)){
  estimates.Comp <- array(dim=c(length(Rs),length(Ds),length(ns),npars,length(qs.comp),length(dmaxs),dimpars))
  times.Comp <- array(dim=c(length(Rs),length(Ds),length(ns),npars,length(qs.comp),length(dmaxs)))
  dimnames(estimates.Comp)<-list(paste('R=',Rs,sep=''),paste('D=',Ds,sep=''),paste('n=',ns,sep=''),paste('par=c(',apply(pars.mat,2,paste,collapse=','),')',sep=''),paste('q.comp=',qs.comp,sep=''),paste('dmax=',dmaxs,sep=''),paste('par ',c(1:dimpars),sep=''))
  dimnames(times.Comp)<-list(paste('R=',Rs,sep=''),paste('D=',Ds,sep=''),paste('n=',ns,sep=''),paste('par=c(',apply(pars.mat,2,paste,collapse=','),')',sep=''),paste('q.comp=',qs.comp,sep=''),paste('dmax=',dmaxs,sep=''))
}
if(any(estimators==3)){
  estimates.Vecchia <- array(dim=c(length(Rs),length(Ds),length(ns),npars,length(qs.vecchia),length(orderings.vecchia),dimpars))
  times.Vecchia <- array(dim=c(length(Rs),length(Ds),length(ns),npars,length(qs.vecchia),length(orderings.vecchia)))
  dimnames(estimates.Vecchia)<-list(paste('R=',Rs,sep=''),paste('D=',Ds,sep=''),paste('n=',ns,sep=''),paste('par=c(',apply(pars.mat,2,paste,collapse=','),')',sep=''),paste('q.vecchia=',qs.vecchia,sep=''),paste('ordering=',orderings.vecchia,sep=''),paste('par ',c(1:dimpars),sep=''))
  dimnames(times.Vecchia)<-list(paste('R=',Rs,sep=''),paste('D=',Ds,sep=''),paste('n=',ns,sep=''),paste('par=c(',apply(pars.mat,2,paste,collapse=','),')',sep=''),paste('q.vecchia=',qs.vecchia,sep=''),paste('ordering=',orderings.vecchia,sep=''))
}

for(r in 1:length(Rs)){
  if(any(estimators==1)){
    estimates.MLE[r,,,,] <- results.list[[r]]$estimates.MLE.r
    times.MLE[r,,,] <- results.list[[r]]$times.MLE.r
  }
  if(any(estimators==2)){
    estimates.Comp[r,,,,,,] <- results.list[[r]]$estimates.Comp.r
    times.Comp[r,,,,,] <- results.list[[r]]$times.Comp.r
  }
  if(any(estimators==3)){
    estimates.Vecchia[r,,,,,,] <- results.list[[r]]$estimates.Vecchia.r
    times.Vecchia[r,,,,,] <- results.list[[r]]$times.Vecchia.r
  }
}

if(any(estimators==1)){
  save(result.list,estimates.MLE,times.MLE,file=paste("data/estimates_MLE_r",paste(range(Rs),collapse="-"),".RData",sep=""))
}
if(any(estimators==2)){
  save(result.list,estimates.Comp,times.Comp,file=paste("data/estimates_Comp_r",paste(range(Rs),collapse="-"),".RData",sep=""))
}
if(any(estimators==3)){
  save(result.list,estimates.Vecchia,times.Vecchia,file=paste("data/estimates_Vecchia_r",paste(range(Rs),collapse="-"),".RData",sep=""))
}





