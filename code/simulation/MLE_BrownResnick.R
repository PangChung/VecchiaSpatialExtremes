#####################################
#####################################
## BROWN-RESNICK MAX-STABLE MODEL ###
#####################################
#####################################

######################
## Full likelihood ###
######################

###########################
## Covariance functions ###
###########################

### Computes covariance matrix
# par: parameter vector
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
Sigma <- function(par,distmat,type=1){
  if(type==1){ ## exponential correlation
    range <- par[1]
    sigma2 <- par[2]
    Sigma.mat <- sigma2*exp(-distmat/range)
    diag(Sigma.mat) <- sigma2
    return( Sigma.mat )
  } else if(type==2){ ## powered exponential correlation
    range <- par[1]
    smooth <- par[2]
    sigma2 <- par[3]
    Sigma.mat <-  sigma2*exp(-(distmat/range)^smooth)
    diag(Sigma.mat) <- sigma2
    return( Sigma.mat )
  } else if(type==3){ ## exchangeable correlation
    corr <- par[1]
    Sigma.mat <- matrix(corr,nrow=nrow(distmat),ncol=ncol(distmat))
    diag(Sigma.mat) <- 1
    return(Sigma.mat)
  }
}

### checks if a given set of parameters yields a valid covariance matrix
# par: parameter vector
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
valid.Sigma <- function(par,distmat,type=1){
  if(type==1){ ## exponential correlation
    range <- par[1]
    sigma2 <- par[2]
    return( (range>0 & sigma2>0) )
  } else if(type==2){ ## powered exponential correlation
    range <- par[1]
    smooth <- par[2]
    sigma2 <- par[3]
    return( (range>0 & smooth>0 & smooth<2 & sigma2>0) )
  } else if(type==3){ ## exchangeable correlation
    corr <- par[1]
    D <- nrow(distmat)
    return( (corr>-1/(D-1) & corr<1) )
  }
}

### returns the lower/upper bounds for all covariance parameters
# par: parameter vector
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
bounds.Sigma <- function(par,distmat,type=1){
  if(type==1){ ## exponential correlation
    range <- par[1]
    sigma2 <- par[2]
    return( rbind(c(0,10*range),c(0,10*sigma2)) )
  } else if(type==2){ ## powered exponential correlation
    range <- par[1]
    smooth <- par[2]
    sigma2 <- par[3]
    return( rbind(c(0,10*range),c(0,2),c(0,10*sigma2)) )
  } else if(type==3){ ## exchangeable correlation
    corr <- par[1]
    D <- nrow(distmat)
    return( matrix(c(-1/(D-1),1),nrow=1,ncol=2) )
  }
}


### Exponent function V for the Brown-Resnick model
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension DxD
V <- function(data,sigma){
   D <- ncol(data)
   if(D==1){
     return(1/data)
   } else{
    fun.i <- function(i){
      eval.i <- t(t(log(data[,-i]/data[,i])) + diag(sigma)[-i]/2 + sigma[i,i]/2 - sigma[i,-i])
      Id <- diag(1,D-1)
      if(i==1){ ### see Wadsworth and Tawn (2014) for the definition of the matrix T.i...
        T.i <- cbind(-1,Id)
      } else if(i==D){
        T.i <- cbind(Id,-1)
      } else{
        T.i <- cbind(Id[,1:(i-1)],-1,Id[,i:(D-1)]) 
      }
      sigma.i <- T.i%*%sigma%*%t(T.i)
      return(apply(eval.i,1,function(x){return(mvtnorm::pmvnorm(upper=x,sigma=sigma.i))})/data[,i])
     }
     return( rowSums(sapply(1:D,fun.i)) )
   }
}

### Negative partial derivatives of the exponent function V for the Brown-Resnick model
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension DxD
# I: vector of indices with respect to which the partial derivatives are computed; if I==c(1:D), the function returns the full mixed derivative.
nVI <- function(data,sigma,I){
  D <- ncol(data)
  nI <- length(I)
  
  if(nI==0){
    return(-V(data,sigma))
  } else if(nI==D){
    if(D==1){
      return(1/data^2)
    } else{
      sigma2 <- sigma[1,1] ## we assume that the process is stationary
    
      sigma.inv <- solve(sigma)
      q <- rowSums(sigma.inv)
      A <- (sigma.inv - q%*%t(q)/sum(q))
    
      log.data <- log(data)
    
      log.Part1 <- 0
      log.Part2 <- ((D-1)/2)*log(2*pi) + (1/2)*log(det(sigma)) + (1/2)*log(sum(q)) + rowSums(log.data)
      log.Part3 <- -(1/2)*(apply(log.data,1,function(x){return(t(x)%*%A%*%x)}) + log.data%*%(2*q/sum(q)) + sigma2 - 1/sum(q))

      res <- exp(log.Part1-log.Part2+log.Part3)
    
      return( drop(res) )
    }
  } else{
    if(D==1){
      return(1/data^2)
    } else{
      sigma2 <- sigma[1,1] ## we assume that the process is stationary
  
      sigma.inv <- solve(sigma)
      q <- rowSums(sigma.inv)
      A <- (sigma.inv - q%*%t(q)/sum(q))
  
      sigma.I <- as.matrix(sigma[I,I])
      sigma.I.inv <- solve(sigma.I)
      q.I <- rowSums(sigma.I.inv)
      A.I <- (sigma.I.inv - q.I%*%t(q.I)/sum(q.I))
  
      K10 <- matrix(0,nrow=D,ncol=nI)
      K10[I,] <- diag(nI)
      K01 <- matrix(0,nrow=D,ncol=D-nI)
      K01[-I,] <- diag(D-nI)
  
      log.data.I <- matrix(log(data[,I]),ncol=nI)
      
      gamma.inv <- t(K01)%*%A%*%K01
      gamma <- chol2inv(chol(gamma.inv))
      mu <- -gamma%*%(t(K01)%*%A%*%K10%*%t(log.data.I) + rowSums(t(K01)%*%sigma.inv)/sum(sigma.inv))
      eval <- log(data[,-I])-t(mu)
      
      log.Part1 <- log( apply(eval,1,function(x){return(mvtnorm::pmvnorm(upper=x,sigma=gamma))}) )
      log.Part2 <- ((nI-1)/2)*log(2*pi) + (1/2)*log(det(sigma.I)) + (1/2)*log(sum(q.I)) + rowSums(log.data.I)
      log.Part3 <- -(1/2)*(apply(log.data.I,1,function(x){return(t(x)%*%A.I%*%x)}) + log.data.I%*%(2*q.I/sum(q.I)) + sigma2 - 1/sum(q.I))
      
      res <- exp(log.Part1-log.Part2+log.Part3)
    
      return( drop(res) )
    }
  }
}

### BIVARIATE Exponent function V for the Brown-Resnick model
# data: matrix of dimension nx2, containing n 2-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension 2x2

V.biv <- function(data,sigma){
  vario <- 2*(sigma[1,1]-sigma[1,2])
  a <- sqrt(vario)
  P12 <- a/2-log(data[,1]/data[,2])/a
  P21 <- a/2-log(data[,2]/data[,1])/a
  return(pnorm(P12)/data[,1]+pnorm(P21)/data[,2])
}

### Negative partial derivatives of the exponent function V for the Brown-Resnick model
# data: matrix of dimension nx2, containing n 2-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension 2x2
# I: vector of indices with respect to which the partial derivatives are computed; if I==c(1:2), the function returns the full mixed derivative.
nVI.biv <- function(data,sigma,I){
  nI <- length(I)
  
  if(nI==0){
    return(-V.biv(data,sigma))
  } else if(nI==2){
    vario <- 2*(sigma[1,1]-sigma[1,2])
    a <- sqrt(vario)
    P12 <- a/2-log(data[,1]/data[,2])/a
    P21 <- a/2-log(data[,2]/data[,1])/a
    return( (dnorm(P12)*(1-P12/a)/data[,1]+dnorm(P21)*(1-P21/a)/data[,2])/(a*data[,1]*data[,2]) )
  } else if(nI==1){
    vario <- 2*(sigma[1,1]-sigma[1,2])
    a <- sqrt(vario)
    P12 <- a/2-log(data[,1]/data[,2])/a
    P21 <- a/2-log(data[,2]/data[,1])/a
    if(I==1){
      return( (pnorm(P12)+dnorm(P12)/a)/data[,1]^2 - dnorm(P21)/(a*data[,1]*data[,2]) )
    } else if(I==2){
      return( (pnorm(P21)+dnorm(P21)/a)/data[,2]^2 - dnorm(P12)/(a*data[,1]*data[,2]) )
    }
  }
}


### Negative log likelihood function for Brown-Resnick data with unit Frechet margins
# par: parameter vector
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
nloglik.BR <- function(par,data,distmat,type=1){
	if(valid.Sigma(par,distmat,type)){
	  #fix random seed (and save the current random seed to restore it at the end)
	  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
	  set.seed(747380)
	  sigma <- Sigma(par,distmat,type)
	  D <- ncol(data)
	  all_combn <- list()
	  for(i in 1:D){
	    all_combn[[i]] <- comb_n(D,i) ## using package `Rfast'
	  }
	  all_nVI <- list() ## will contain all the terms nVI (total number is equal to 2^D-1), used later to assemble the log-likelihood...
	  for(i in 1:D){
	    all_nVI[[i]] <- matrix(nrow=nrow(data),ncol=ncol(all_combn[[i]]))
	    for(j in 1:ncol(all_combn[[i]])){
	      all_nVI[[i]][,j] <- nVI(data,sigma,all_combn[[i]][,j])
	    }
	  }
	  get.nVI <- function(I){
	    nI <- length(I)
	    return(all_nVI[[nI]][,which(apply(all_combn[[nI]],2,function(x){return(all(I%in%x))}))])
	  }
	  parts <- listParts(length(c(1:D))) ## using package `partitions'
	  contribution.partition <- function(partition){
	    return( rowprods(as.matrix(as.data.frame(lapply(partition,FUN=function(I){return(get.nVI(I))})))) )
	  }
		res <- sum( log(rowSums(as.matrix(as.data.frame(lapply(parts,contribution.partition))))) - V(data,sigma) )
		
		#restore random seed to its previous value
		assign(".Random.seed", oldSeed, envir=globalenv())
		
		return(-res)
	} else{
		return(Inf)	
	}
}

### Negative log likelihood function for BIVARIATE Brown-Resnick data with unit Frechet margins
# par: parameter vector
# data: matrix of dimension nx2, containing n 2-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: matrix of dimension 2x2, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
nloglik.BR.biv <- function(par,data,distmat,type=1){
  if(valid.Sigma(par,distmat,type)){
    #fix random seed (and save the current random seed to restore it at the end)
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    
    sigma <- Sigma(par,distmat,type)
    
    res <- sum( log(nVI.biv(data,sigma,1)*nVI.biv(data,sigma,2)+nVI.biv(data,sigma,1:2))-V.biv(data,sigma) )
    
    #restore random seed to its previous value
    assign(".Random.seed", oldSeed, envir=globalenv())
    
    return(-res)
  } else{
    return(Inf)	
  }
}


### Function that returns the MLE for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
MLE.BR <- function(data,init,fixed,distmat,type){
  t <- proc.time()
  nloglik.BR2 <- function(par2,data,distmat,type){
    par <- init
    par[!fixed] <- par2
    return(nloglik.BR(par,data,distmat,type))
  }
  init2 <- init[!fixed]
  if(length(init2)==1){
    bounds.Sig <- bounds.Sigma(init,distmat,type)[!fixed,]
    opt <- optim(par=init2,fn=nloglik.BR2,data=data,distmat=distmat,type=type,method="Brent",lower=bounds.Sig[1],upper=bounds.Sig[2],control=list(maxit=10000))
  } else{
    opt <- optim(par=init2,fn=nloglik.BR2,data=data,distmat=distmat,type=type,method="Nelder-Mead",control=list(maxit=10000))
  }
  time <- proc.time()-t
  res <- opt
  res$time <- time[3]
	return( res )
}

### Function that returns the MLE for the BIVARIATE Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nx2, containing n 2-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: matrix of dimension 2x2, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
MLE.BR.biv <- function(data,init,fixed,distmat,type){
  t <- proc.time()
  nloglik.BR.biv2 <- function(par2,data,distmat,type){
    par <- init
    par[!fixed] <- par2
    return(nloglik.BR.biv(par,data,distmat,type))
  }
  init2 <- init[!fixed]
  if(length(init2)==1){
    bounds.Sig <- bounds.Sigma(init,distmat,type)[!fixed,]
    opt <- optim(par=init2,fn=nloglik.BR.biv2,data=data,distmat=distmat,type=type,method="Brent",lower=bounds.Sig[1],upper=bounds.Sig[2],control=list(maxit=10000))
  } else{
    opt <- optim(par=init2,fn=nloglik.BR.biv2,data=data,distmat=distmat,type=type,method="Nelder-Mead",control=list(maxit=10000))
  }
  time <- proc.time()-t
  res <- opt
  res$time <- time[3]
  return( res )
}



###########################
## Composite likelihood ###
###########################

### Negative log composite-likelihood function for Brown-Resnick data with unit Frechet margins
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# index: q-by-Q matrix of q-dimensional margins to be used in the composite likelihood. Here Q refers to the number of composite likelihood contributions (with 1<=Q<=choose(D,q))
nlogcomplik.BR <- function(par,data,distmat,type,index){
  if(valid.Sigma(par,distmat,type)){
    nlogcomplik.contribution.BR <- function(index){
      return(nloglik.BR(par,data[,index],distmat[index,index],type))
    }
    res <- sum(apply(index,2,nlogcomplik.contribution.BR))
    return(res)
  } else{
    return(Inf)	
  }
}

### Negative log pairwise-likelihood function for Brown-Resnick data with unit Frechet margins
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# index: 2-by-Q matrix of 2-dimensional margins to be used in the pairwise likelihood. Here Q refers to the number of pairwise likelihood contributions (with 1<=Q<=choose(D,2))
nlogcomplik.BR.biv <- function(par,data,distmat,type,index){
  if(valid.Sigma(par,distmat,type)){
    nlogcomplik.contribution.BR <- function(index){
      return(nloglik.BR.biv(par,data[,index],distmat[index,index],type))
    }
    res <- sum(apply(index,2,nlogcomplik.contribution.BR))
    return(res)
  } else{
    return(Inf)	
  }
}

### Function that returns the MCLE (maximum composite likelihood estimator) for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# index: q-by-Q matrix of q-dimensional margins to be used in the composite likelihood. Here Q refers to the number of composite likelihood contributions (with 1<=Q<=choose(D,q)).
MCLE.BR <- function(data,init,fixed,distmat,type,index){
  t <- proc.time()
  nlogcomplik.BR2 <- function(par2,data,distmat,type,index){
    par <- init
    par[!fixed] <- par2
    return(nlogcomplik.BR(par,data,distmat,type,index))
  }
  init2 <- init[!fixed]
  if(length(init2)==1){
    bounds.Sig <- bounds.Sigma(init,distmat,type)[!fixed,]
    opt <- optim(par=init2,fn=nlogcomplik.BR2,data=data,distmat=distmat,type=type,index=index,method="Brent",lower=bounds.Sig[1],upper=bounds.Sig[2],control=list(maxit=10000))
  } else{
    opt <- optim(par=init2,fn=nlogcomplik.BR2,data=data,distmat=distmat,type=type,index=index,method="Nelder-Mead",control=list(maxit=10000))
  }
  time <- proc.time()-t
  res <- opt
  res$time <- time[3]
  return( res )
}

### Function that returns the MCLE (maximum pairwise likelihood estimator) for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: matrix of dimension 2x2, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# index: 2-by-Q matrix of 2-dimensional margins to be used in the pairwise likelihood. Here Q refers to the number of composite likelihood contributions (with 1<=Q<=choose(D,2)).
MCLE.BR.biv <- function(data,init,fixed,distmat,type,index){
  t <- proc.time()
  nlogcomplik.BR.biv2 <- function(par2,data,distmat,type,index){
    par <- init
    par[!fixed] <- par2
    return(nlogcomplik.BR.biv(par,data,distmat,type,index))
  }
  init2 <- init[!fixed]
  if(length(init2)==1){
    bounds.Sig <- bounds.Sigma(init,distmat,type)[!fixed,]
    opt <- optim(par=init2,fn=nlogcomplik.BR.biv2,data=data,distmat=distmat,type=type,index=index,method="Brent",lower=bounds.Sig[1],upper=bounds.Sig[2],control=list(maxit=10000))
  } else{
    opt <- optim(par=init2,fn=nlogcomplik.BR.biv2,data=data,distmat=distmat,type=type,index=index,method="Nelder-Mead",control=list(maxit=10000))
  }
  time <- proc.time()-t
  res <- opt
  res$time <- time[3]
  return( res )
}



##############################
## Vecchia's approximation ###
##############################

### Negative log likelihood function for Brown-Resnick data with unit Frechet margins, based on Vecchia's approximation
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an q-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where q is the number of neighbours, i.e., the size of the conditioning set)
nlogVecchialik.BR <- function(par,data,distmat,type,vecchia.seq,neighbours){
  if(valid.Sigma(par,distmat,type)){
    logVecchialik.contribution.BR <- function(i){
      if(i==1){
        contribution <- nloglik.BR(par,as.matrix(data[,vecchia.seq[1]]),distmat=as.matrix(distmat[vecchia.seq[1],vecchia.seq[1]]),type) #density of 1st variable in the sequence (unit Fréchet)
      } else{
        ind.i <- vecchia.seq[i] #index of ith-variable in the Vecchia sequence
        ind.neighbours <- na.omit(neighbours[,i])
        num <- nloglik.BR(par,as.matrix(data[,c(ind.i,ind.neighbours)]),as.matrix(distmat[c(ind.i,ind.neighbours),c(ind.i,ind.neighbours)]),type) #joint density of ith-variable and its conditioning set
        denom <- nloglik.BR(par,as.matrix(data[,c(ind.neighbours)]),as.matrix(distmat[c(ind.neighbours),c(ind.neighbours)]),type) #joint density of conditioning set only
        contribution <- num-denom
      }
      return(contribution)
    }
    res <- sum(sapply(c(1:length(vecchia.seq)),FUN=logVecchialik.contribution.BR))
    return(res)
  } else{
    return(Inf)	
  }
}

### Negative log likelihood function for Brown-Resnick data with unit Frechet margins, based on Vecchia's approximation based on PAIRWISE margins
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an 1-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where 1 is the number of neighbours, i.e., the size of the conditioning set)
nlogVecchialik.BR.biv <- function(par,data,distmat,type,vecchia.seq,neighbours){
  if(valid.Sigma(par,distmat,type)){
    logVecchialik.contribution.BR <- function(i){
      if(i==1){
        contribution <- nloglik.BR(par,as.matrix(data[,vecchia.seq[1]]),distmat=as.matrix(distmat[vecchia.seq[1],vecchia.seq[1]]),type) #density of 1st variable in the sequence (unit Fréchet)
      } else{
        ind.i <- vecchia.seq[i] #index of ith-variable in the Vecchia sequence
        ind.neighbours <- na.omit(neighbours[,i])
        num <- nloglik.BR.biv(par,as.matrix(data[,c(ind.i,ind.neighbours)]),as.matrix(distmat[c(ind.i,ind.neighbours),c(ind.i,ind.neighbours)]),type) #joint density of ith-variable and its conditioning set
        denom <- nloglik.BR(par,as.matrix(data[,c(ind.neighbours)]),as.matrix(distmat[c(ind.neighbours),c(ind.neighbours)]),type) #joint density of conditioning set only
        contribution <- num-denom
      }
      return(contribution)
    }
    res <- sum(sapply(c(1:length(vecchia.seq)),FUN=logVecchialik.contribution.BR))
    return(res)
  } else{
    return(Inf)	
  }
}

### Function that returns the MVLE (maximum Vecchia likelihood estimator) for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an q-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where q is the number of neighbours, i.e., the size of the conditioning set)
MVLE.BR <- function(data,init,fixed,distmat,type,vecchia.seq,neighbours){
  t <- proc.time()
  nlogVecchialik.BR2 <- function(par2,data,distmat,type,vecchia.seq,neighbours){
    par <- init
    par[!fixed] <- par2
    return(nlogVecchialik.BR(par,data,distmat,type,vecchia.seq,neighbours))
  }
  init2 <- init[!fixed]
  if(length(init2)==1){
    bounds.Sig <- bounds.Sigma(init,distmat,type)[!fixed,]
    opt <- optim(par=init2,fn=nlogVecchialik.BR2,data=data,distmat=distmat,type=type,vecchia.seq=vecchia.seq,neighbours=neighbours,method="Brent",lower=bounds.Sig[1],upper=bounds.Sig[2],control=list(maxit=10000))
  } else{
    opt <- optim(par=init2,fn=nlogVecchialik.BR2,data=data,distmat=distmat,type=type,vecchia.seq=vecchia.seq,neighbours=neighbours,method="Nelder-Mead",control=list(maxit=10000))
  }
  time <- proc.time()-t
  res <- opt
  res$time <- time[3]
  return( res )
}

### Function that returns the MVLE (maximum Vecchia likelihood estimator) for the Brown-Resnick model with unit Frechet margins, using PAIRWISE components
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: matrix of dimension DxD, containing all pairwise distances between locations
# type: choose between different types of covariance matrices
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an 1-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where 1 is the number of neighbours, i.e., the size of the conditioning set)
MVLE.BR.biv <- function(data,init,fixed,distmat,type,vecchia.seq,neighbours){
  t <- proc.time()
  nlogVecchialik.BR.biv2 <- function(par2,data,distmat,type,vecchia.seq,neighbours){
    par <- init
    par[!fixed] <- par2
    return(nlogVecchialik.BR.biv(par,data,distmat,type,vecchia.seq,neighbours))
  }
  init2 <- init[!fixed]
  if(length(init2)==1){
    bounds.Sig <- bounds.Sigma(init,distmat,type)[!fixed,]
    opt <- optim(par=init2,fn=nlogVecchialik.BR.biv2,data=data,distmat=distmat,type=type,vecchia.seq=vecchia.seq,neighbours=neighbours,method="Brent",lower=bounds.Sig[1],upper=bounds.Sig[2],control=list(maxit=10000))
  } else{
    opt <- optim(par=init2,fn=nlogVecchialik.BR.biv2,data=data,distmat=distmat,type=type,vecchia.seq=vecchia.seq,neighbours=neighbours,method="Nelder-Mead",control=list(maxit=10000))
  }
  time <- proc.time()-t
  res <- opt
  res$time <- time[3]
  return( res )
}



