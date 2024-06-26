###########################################################
###########################################################
## BROWN-RESNICK MAX-STABLE MODEL with General Variogram###
###########################################################
###########################################################

######################
## Full likelihood ###
######################

###########################
## Variogram functions ####
###########################
##input: values in R and need re-parametrization
vario.func <- function(loc,par){ ##return a covariance matrix
  alpha = par[1];lambda = par[2];a = par[3]; theta = par[4]
  #Sigma <- matrix(c(par[3],-par[4],-par[4],1),2,2)
  A = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
  Sigma <- A%*%diag(c(1,a),2)%*%t(A)
  loc = matrix(loc,ncol=2)
  n = nrow(loc)
  if(n==1){
    val=2*(sqrt(t(loc[1,])%*%Sigma%*%loc[1,])/lambda)^alpha
    return(val)
  }
  fun <- function(idx){
    loc.temp <- loc[idx,]
      if(idx[1]==idx[2]){
        h <- loc.temp[1,]
        val = 2*(sqrt(t(h)%*%Sigma%*%h)/lambda)^alpha
      }else{
        h <- loc.temp[1,]-loc.temp[2,]
        val <- (sqrt(t(loc.temp[1,])%*%Sigma%*%(loc.temp[1,]))/lambda)^alpha+
          (sqrt(t(loc.temp[2,])%*%Sigma%*%(loc.temp[2,]))/lambda)^alpha-
          (sqrt(t(h)%*%Sigma%*%h)/lambda)^alpha 
      }
      return(val)
    }
  idx <- cbind(rep(1:n, times = n:1),unlist(lapply(1:n, function(x){x:n})))
  val <- apply(idx,1,fun)
  val.mat <- matrix(1,n,n)
  val.mat[idx]<-val;
  val.mat[idx[,c(2,1)]]<-val
  return(val.mat + .Machine$double.eps * diag(n))
}

# check the paramaters 
par.check <- function(par){
  return( (par[1] > 0 & par[1] < 2 & par[2] > 0.01 & par[2] < 1000) )
}

# par.check <- function(par){
#   return( (par[1] > 0 & par[2] > 10 & par[2] < 1000 & par[3] > 0 & par[4] >= 0 ))
# }


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
  if(nI==0){ ## If the index set is empty
    return(-V(data,sigma))
  }else if(nI==D){ ## full derivative
    if(D==1){
      return(1/data^2)
    } else{
      sigma.DD <- diag(sigma)
      sigma.inv <- chol2inv(chol(sigma))
      q <- rowSums(sigma.inv)
      q.sum <- sum(q)
      A <- (sigma.inv - q%*%t(q)/q.sum)
      sigma.q <- t(sigma.DD)%*%q
      
      log.data <- log(data)
      
      log.Part1 <- 0
      log.Part2 <- ((D-1)/2)*log(2*pi) + (1/2)*log(det(sigma)) + (1/2)*log(q.sum) + rowSums(log.data)
      log.Part3 <- c(-(1/2)*( 1/4*t(sigma.DD)%*%sigma.inv%*%sigma.DD -1/4*(sigma.q)^2/q.sum + sigma.q/q.sum - 1/q.sum))
      log.Part4 <- c(-(1/2)*(apply(log.data,1,function(x){return(t(x)%*%A%*%x)}) + log.data%*%(q%*%(2-sigma.q)/q.sum + sigma.inv%*%sigma.DD)))
      
      res <- exp(log.Part1-log.Part2+log.Part3+log.Part4)
      return( drop(res) )
    }
  } else{ ## Partial derivative
    if(D==1){ ## If there is only one variable
      return(1/data^2)
    } else{
        sigma.DD <- diag(sigma)
        sigma.II <- sigma.DD[I]
        sigma.inv <- chol2inv(chol(sigma))
        q <- rowSums(sigma.inv)
        A <- (sigma.inv - q%*%t(q)/sum(q))
        
        
        sigma.I <- as.matrix(sigma[I,I])
        sigma.I.inv <- chol2inv(chol(sigma.I))
        q.I <- rowSums(sigma.I.inv)
        A.I <- (sigma.I.inv - q.I%*%t(q.I)/sum(q.I))
        
        K10 <- matrix(0,nrow=D,ncol=nI)
        K10[I,] <- diag(nI)
        K01 <- matrix(0,nrow=D,ncol=D-nI)
        K01[-I,] <- diag(D-nI)
        
        log.data <- log(data)
        log.data.I <- matrix(log.data[,I],ncol=nI)
        
        gamma.inv <- t(K01)%*%A%*%K01
        gamma <- chol2inv(chol(gamma.inv))
        mu <- -gamma%*%(t(K01)%*%A%*%K10%*%t(log.data.I) + c(t(K01)%*%( (q - 1/2*q%*%t(q)%*%sigma.DD )/sum(q) + 1/2*sigma.inv%*%sigma.DD )) )
        eval.x <- log(data[,-I])-t(mu)
        sigma.q.II <- t(sigma.II)%*%q.I
        q.I.sum <- sum(q.I)
        
        log.Part1 <- apply(eval.x,1,function(x){return(max(pmvnorm(upper=x,sigma=gamma),0))})
        log.Part2 <- ((nI-1)/2)*log(2*pi) + (1/2)*log(det(sigma.I)) + (1/2)*log(q.I.sum) + rowSums(log.data.I)
        log.Part3 <- c(-(1/2)*( 1/4*t(sigma.II)%*%sigma.I.inv%*%sigma.II -1/4*(sigma.q.II)^2/q.I.sum + sigma.q.II/q.I.sum - 1/q.I.sum))
        log.Part4 <- c(-(1/2)*(apply(log.data.I,1,function(x){return(t(x)%*%A.I%*%x)}) + log.data.I%*%(q.I%*%(2-sigma.q.II)/q.I.sum + sigma.I.inv%*%sigma.II)))
        res <- drop(log.Part1*exp(-log.Part2+log.Part3+log.Part4))
      return(drop(res))
    }
  }
}

### BIVARIATE Exponent function V for the Brown-Resnick model
# data: matrix of dimension nx2, containing n 2-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension 2x2

V.biv <- function(data,sigma){
  vario <- (sigma[1,1]+sigma[2,2]-2*sigma[1,2])
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
    vario <- (sigma[1,1]+sigma[2,2]-2*sigma[1,2])
    a <- sqrt(vario)
    P12 <- a/2-log(data[,1]/data[,2])/a
    P21 <- a/2-log(data[,2]/data[,1])/a
    return( (dnorm(P12)*(1-P12/a)/data[,1]+dnorm(P21)*(1-P21/a)/data[,2])/(a*data[,1]*data[,2]) )
  } else if(nI==1){
    vario <- (sigma[1,1]+sigma[2,2]-2*sigma[1,2])
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
# distmat: matrix of dimension DxD, containing all pairwise distances between locations or Locations
# FUN: the variogram function that returns the covraiance matrix
nloglik.BR <- function(par,data,distmat,FUN){
    #fix random seed (and save the current random seed to restore it at the end)
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma <- FUN(distmat,par)
    # sigma = FUN(par,distmat)
    D <- ncol(data)
    all_combn <- lapply(1:D,FUN=combn,x=D,simplify=FALSE) ## not using package `Rfast' and return a list of lists
    all_nVI <- list() ## will contain all the terms nVI (total number is equal to 2^D-1), used later to assemble the log-likelihood...
    tryCatch({
    all_nVI <- lapply(all_combn,FUN = function(idx){vapply(idx,nVI,FUN.VALUE = rep(0,nrow(data)),data=data,sigma=sigma)})
    get.nVI <- function(I){
      nI <- length(I)
      return(all_nVI[[nI]][,which(sapply(all_combn[[nI]],function(x){return(all(I%in%x))}))])
    }
    parts <- listParts(length(c(1:D))) ## using package `partitions'
    contribution.partition <- function(partition){
      return( apply(as.matrix(as.data.frame(lapply(partition,FUN=get.nVI))),1,prod) )
    }
    res <- log(rowSums(as.matrix(as.data.frame(lapply(parts,contribution.partition))))) - V(data,sigma)
    res <- mean(res[is.finite(res)])
    #restore random seed to its previous value
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(-res)},error=function(e){NA})
}

###########################
## Composite likelihood ###
###########################

### Negative log composite-likelihood function for Brown-Resnick data with unit Frechet margins
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: coordinates
# FUN: function returns covariance matrix
# index: q-by-Q matrix of q-dimensional margins to be used in the composite likelihood. Here Q refers to the number of composite likelihood contributions (with 1<=Q<=choose(D,q))
nlogcomplik.BR <- function(par,data,distmat,FUN,index,ncores){
    nlogcomplik.contribution.BR <- function(index){
      val <- nloglik.BR(par,data[,index],distmat[index,],FUN)
    }
    res <- mean(unlist(mclapply(as.list(as.data.frame(index)),nlogcomplik.contribution.BR,mc.cores = ncores,mc.set.seed = F)),na.rm = TRUE)
    return(res)
}


### Function that returns the MCLE (maximum composite likelihood estimator) for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: coordinates
# FUN: function returns covariance matrix
# index: q-by-Q matrix of q-dimensional margins to be used in the composite likelihood. Here Q refers to the number of composite likelihood contributions (with 1<=Q<=choose(D,q)).
MCLE.BR <- function(data,init,fixed,distmat,FUN,index,ncores,maxit=100,method="Nelder-Mead",hessian=FALSE){
  t <- proc.time()
  nlogcomplik.BR2 <- function(par2){
    par <- init2
    par[!fixed] <- par2
    if(!par.check(par)){return(Inf)}
    return(nlogcomplik.BR(par,data,distmat,FUN,index,ncores))
  }
  fixed0 <- fixed
  index0 <- index
  init2 <- init
  lower = c(0.1,0.1,0,-Inf)
  upper = c(1.9,1000,Inf,Inf)
  val_fn = c()
  if(sum(!fixed)!=1){
    opt <- optim(par=init2[!fixed],fn=nlogcomplik.BR2,method="Nelder-Mead",control=list(maxit=maxit,trace=TRUE),hessian=hessian)
    init2[!fixed] = opt$par
  }
  if(sum(!fixed)==1){
    opt <- optim(par=init2[!fixed],fn=nlogcomplik.BR2,lower = lower[!fixed],upper=upper[!fixed], method="Brent",control=list(maxit=maxit,trace=TRUE),hessian=hessian)
    init2[!fixed] <- opt$par 
  }
  time <- proc.time()-t
  res <- opt
  res$par = init2
  res$time <- time[3]
  return( res )
}

##############################
## Vecchia's approximation ###
##############################

### Negative log likelihood function for Brown-Resnick data with unit Frechet margins, based on Vecchia's approximation
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# distmat: coordinates
# FUN: function returns covariance matrix
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an q-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where q is the number of neighbours, i.e., the size of the conditioning set)
nlogVecchialik.BR <- function(par,data,distmat,FUN,vecchia.seq,neighbours,ncores){
    logVecchialik.contribution.BR <- function(i){
      if(i==1 & !any(!is.na(neighbours[,i]))){
        contribution <- nloglik.BR(par,as.matrix(data[,vecchia.seq[1]]),distmat=matrix(distmat[vecchia.seq[1],],ncol=2),FUN) #density of 1st variable in the sequence (unit Fréchet)
      } else{
        ind.i <- vecchia.seq[i] #index of ith-variable in the Vecchia sequence
        ind.neighbours <- na.omit(neighbours[,i])
        num <- nloglik.BR(par,as.matrix(data[,c(ind.i,ind.neighbours)]),matrix(distmat[c(ind.i,ind.neighbours),],ncol=2),FUN) #joint density of ith-variable and its conditioning set
        denom <- nloglik.BR(par,as.matrix(data[,c(ind.neighbours)]),matrix(distmat[c(ind.neighbours),],ncol=2),FUN) #joint density of conditioning set only
        contribution <- num-denom
      }
      return(contribution)
    }
    res <- mean(unlist(mclapply(1:length(vecchia.seq),FUN=logVecchialik.contribution.BR,mc.cores=ncores,mc.set.seed = F)),na.rm=TRUE)
    return(res)
}

### Function that returns the MVLE (maximum Vecchia likelihood estimator) for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# distmat: coordinates
# FUN: function returns covariance matrix
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an q-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where q is the number of neighbours, i.e., the size of the conditioning set)
MVLE.BR <- function(data,init,fixed,distmat,FUN,vecchia.seq,neighbours,ncores,maxit=100,method="Nelder-Mead",hessian=FALSE){
  t <- proc.time()
  nlogVecchialik.BR2 <- function(par2){
    par <- init2
    par[!fixed] <- par2
    if(!par.check(par)){return(Inf)}
    return(nlogVecchialik.BR(par,data,distmat,FUN,vecchia.seq,neighbours,ncores))
  }
  
  fixed0 <- fixed
  init2 <- init
  lower = c(0.1,0.1,0,-Inf)
  upper = c(1.5,1000,Inf,Inf)
  val_fn = c()
  if(sum(!fixed)!=1){
    opt <- optim(par=init2[!fixed],fn=nlogVecchialik.BR2,method="Nelder-Mead",control=list(maxit=maxit,trace=TRUE),hessian=hessian)
    init2[!fixed] = opt$par
  }
  if(sum(!fixed)==1){
      opt <- optim(par=init2[!fixed],fn=nlogVecchialik.BR2,lower = lower[!fixed],upper=upper[!fixed],method="Brent",control=list(maxit=maxit,trace=TRUE),hessian=hessian)
      init2[!fixed] <- opt$par
  }
  time <- proc.time()-t
  res = opt
  res$par = init2
  res$time <- time[3]
  return( res )
}

neighbours <- function(ind,vecchia.seq,q,distmat){
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

FitVecchia <- function(data,loc,init,fixed,vecchia.seq,q,vario=vario.func){
	D = nrow(loc)

	neighbours.mat <- sapply(1:D,FUN=neighbours,vecchia.seq=vecchia.seq,q=q)

	if(q==1){ neighbours.mat <- matrix(neighbours.mat,nrow=1) }

	fit.result <- MVLE.BR(data=data,init=init,fixed=fixed,
	distmat=loc,FUN = vario,vecchia.seq=vecchia.seq,
	neighbours=neighbours.mat,ncores,method="Nelder-Mead",maxit=1000,hessian=hessian)

	return(fit.result)
}

