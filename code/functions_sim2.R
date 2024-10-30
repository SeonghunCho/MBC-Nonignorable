########################################################
## Functions
########################################################
Estphi_MAR <- function(Zmat,deltavec){
  phihat <- rep(0,ncol(Zmat))
  count <- 0
  while(TRUE){
    count <- count+1
    phiprev <- phihat
    etavec <- c(Zmat%*%phihat)
    pivec <- 1/(1+exp(-etavec))
    ll <- sum(deltavec*log(pivec) + (1-deltavec)*log(1-pivec))
    dll <- c(t(Zmat)%*%(deltavec-pivec))
    d2ll <- - t(Zmat)%*%diag( pivec*(1-pivec) )%*%Zmat
    phiup <- -c(ginv(d2ll,tol=.Machine$double.eps)%*%dll)
    d <- sqrt(sum(phiup^2))
    # cat(count," : d = ",d,", ll = ",ll,"\n",sep="")
    if(d < 1e-6) break
    step_size <- 1
    while(TRUE){
      phinew <- phihat+step_size*phiup
      etavec <- c(Zmat%*%phinew)
      pivec <- 1/(1+exp(-etavec))
      if(any(pivec==0)|any(pivec==1)){
        step_size <- 0.5*step_size
        next
      }
      ll_new <- sum(deltavec*log(pivec) + (1-deltavec)*log(1-pivec))
      if(ll_new<ll){
        step_size <- 0.5*step_size
      }else{
        break
      }
    }
    phihat <- phinew
  }
  return(phihat)
}
Estphi <- function(Bmat,Zsmat,deltavec){
  n <- length(deltavec)
  n1 <- sum(deltavec)
  n0 <- n-n1
  Bsmat <- Bmat[deltavec,]
  phihat <- rep(0,ncol(Zsmat))
  q <- length(phihat)
  count <- 0
  while(TRUE){
    count <- count+1
    phiprev <- phihat
    pisvec <- 1/(1+exp(-c(Zsmat%*%phihat)))
    pivec <- rep(1,n)
    pivec[deltavec] <- pisvec
    U <- t(Bmat)%*%( deltavec/pivec-1 )
    dU <- -t(Bsmat)%*%diag( (1-pisvec)/pisvec )%*%Zsmat
    phiup <- -c(ginv(t(dU)%*%dU,tol=.Machine$double.eps)%*%t(dU)%*%U)
    d <- sqrt(sum(phiup^2))
    # cat(count," : d = ",d,"\n",sep="")
    if(d < 1e-6) break
    step_size <- 1
    while(TRUE){
      phinew <- phihat+step_size*phiup
      pisvec <- 1/(1+exp(-c(Zsmat%*%phinew)))
      if(any(pisvec==0)){
        step_size <- step_size*0.5
        next
      }
      pivec <- rep(1,n)
      pivec[deltavec] <- pisvec
      Unew <- t(Bmat)%*%( deltavec/pivec-1 )
      if(sum(Unew^2)>=sum(U^2)){
        step_size <- step_size*0.5
      }else{
        break
      }
    }
    phihat <- phinew
    # if(count>=100){
    #   Suc <- FALSE
    #   break
    # }
  }
  return(phihat)
}
Estquan <- function(quan,ysvec,pvec,alpha){
  return( sum(pvec*(ysvec<=quan))/sum(pvec)-alpha )
}
########################################################
#    This part is originated from Owen.
#    There are some differences between this version and Owen's original version.
#
#    The function llog() is equal to the natural
#  logarithm on the interval from eps >0 to infinity.
#  Between -infinity and eps, llog() is a quadratic.
#  llogp() and llogpp() are the first two derivatives
#  of llog().  All three functions are continuous
#  across the "knot" at eps.
#
#    A variation with a second knot at a large value
#  M did not appear to work as well.
#
#    The cutoff point, eps, is usually 1/n, where n
#  is the number of observations.  Unless n is extraordinarily
#  large, dividing by eps is not expected to cause numerical
#  difficulty.
llog <- function( z, eps ){
  ans <- z
  lo <- (z<eps)
  ans[ lo  ] <- log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
  ans[ !lo ] <- log( z[!lo] )
  return(ans)
}
llogp <- function( z, eps ){
  ans <- z
  lo <- (z<eps)
  ans[ lo  ] <- 2.0/eps - z[lo]/eps^2
  ans[ !lo ] <- 1/z[!lo]
  return(ans)
}
llogpp <- function( z, eps ){
  ans <- z
  lo <- (z<eps)
  ans[ lo  ] <- -1.0/eps^2
  ans[ !lo ] <- -1.0/z[!lo]^2
  return(ans)
}
logelr <- function( P, A, n0, lam, eps){
  P <- as.matrix(P)
  A <- as.matrix(A)
  n1 <- nrow(P)
  K <- ncol(P)
  n <- n1+n0
  z <- cbind(P-1,A)
  arg <- as.vector( n/n1 + z %*% lam )
  return( -sum( llog(arg,eps) ) - n0*llog(sum(lam[1:K]),eps) )
}
EstEL <- function(P,A,n0,eps,maxit=200){
  P <- as.matrix(P)
  A <- as.matrix(A)
  n1 <- nrow(P)
  K <- ncol(P)
  r <- ncol(A)
  
  n <- n1+n0
  z <- cbind(P-1,A)
  
  TINY <- sqrt( .Machine$double.xmin )
  svdtol <- TINY
  gradtol <- TINY
  
  lam <- c( rep(n/(K*n1),K) , rep(0,r) )
  
  #
  #    Preset the weights for combining Newton and gradient
  # steps at each of 16 inner iterations, starting with
  # the Newton step and progressing towards shorter vectors
  # in the gradient direction.  Most commonly only the Newton
  # step is actually taken, though occasional step reductions
  # do occur.
  #
  
  nwts <- c( 3^-c(0:3), rep(0,12) )
  gwts <- 2^( -c(0:(length(nwts)-1)))
  gwts <- (gwts^2 - nwts^2)^.5
  gwts[12:16] <- gwts[12:16] * 10^-c(1:5)
  
  #
  #    Iterate, finding the Newton and gradient steps, and
  # choosing a step that reduces the objective if possible.
  #
  
  nits <- 0
  gsize <- gradtol + 1
  while(  nits<maxit && gsize > gradtol  ){
    sl <- sum(lam[1:K])
    arg <- as.vector( n/n1 + z %*% lam )
    wts1 <- as.vector( llogp(arg, eps) )
    wts2 <- as.vector( -llogpp(arg, eps) )
    grad <- -t(z)%*%wts1
    grad[1:K] <- grad[1:K] - n0*llogp(sl,eps)
    gsize <- mean( abs(grad) )
    
    Hess <- t(z*wts2)%*%z
    Hess[1:K,1:K] <- Hess[1:K,1:K] + n0*llogpp(sl,eps)
    eigenH <- eigen( Hess )
    if( min(eigenH$values) < max(eigenH$values)*svdtol )
      eigenH$values <- eigenH$values + max(eigenH$values)*svdtol
    nstep <- eigenH$vectors %*% (t(eigenH$vectors)/eigenH$values)
    nstep <- as.vector( nstep %*% grad )
    
    gstep <- -grad
    if(  sum(nstep^2) < sum(gstep^2) )
      gstep <- gstep*sum(nstep^2)^.5/sum(gstep^2)^.5
    ologelr <- logelr(P, A, n0, lam, eps)
    ninner <- 0
    for(  i in 1:length(nwts) ){
      nlogelr <- logelr( P, A, n0, lam+nwts[i]*nstep+gwts[i]*gstep, eps )
      if( nlogelr < ologelr ){
        lam <- lam+nwts[i]*nstep+gwts[i]*gstep
        ninner <- i
        break
      }
    }
    nits <- nits+1
    if( ninner==0 ) nits <- maxit
  }
  arg <- as.vector( n/n1 + z %*% lam )
  pvec <- as.vector( llogp(arg, eps) )
  pvec <- pvec/sum(pvec)
  return(pvec)
}
########################################################
CCA <- function(Xmat,yvec,deltavec,alpha){
  Xsmat <- Xmat[deltavec,]
  ysvec <- yvec[deltavec]
  betahat <- c(solve(t(Xsmat)%*%Xsmat)%*%t(Xsmat)%*%ysvec)
  muhat <- mean(ysvec)
  quanhat <- quantile(ysvec,prob=alpha,names=F)
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
VAL <- function(Xval,yval,alpha){
  nval <- length(yval)
  betahat <- c(solve(t(Xval)%*%Xval)%*%t(Xval)%*%yval)
  muhat <- mean(yval)
  quanhat <- quantile(yval,prob=alpha,names=F)
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
VALnuis <- function(Xval,yval){
  nval <- length(yval)
  betahat <- c(solve(t(Xval)%*%Xval)%*%t(Xval)%*%yval)
  nuishat <- c(betahat)
  return(nuishat)
}
CEL <- function(Xmat,yvec,deltavec,Aux,is.mar,Zmat,Zsmat,Bmat,alpha){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  Xsmat <- Xmat[deltavec,]
  ysvec <- yvec[deltavec]
  ##########################################
  if(is.mar){
    phihat <- Estphi_MAR(Zmat,deltavec)
  }else{
    phihat <- Estphi(Bmat,Zsmat,deltavec)
  }
  pivec <- 1/(1+exp(-c(Zsmat%*%phihat)))
  ##########################################
  pimat <- matrix(pivec,ncol=1)
  pvec <- EstEL(pimat,Aux,n-n1,1/(10*n))
  betahat <- c(solve(t(Xsmat)%*%diag(pvec)%*%Xsmat)%*%t(Xsmat)%*%diag(pvec)%*%ysvec)
  muhat <- sum(pvec*ysvec)/sum(pvec)
  res_quan <- uniroot(Estquan,interval=range(ysvec),ysvec=ysvec,pvec=pvec,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
MCEL <- function(nuis,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,Bmat,alpha){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  beta <- nuis
  s2 <- s20
  yhat <- c(Xmat%*%beta)
  Xsmat <- Xmat[deltavec,]
  ysvec <- yvec[deltavec]
  ##########################################
  K <- length(is.mar)
  pimat <- NULL
  for(k in 1:K){
    if(is.mar[k]){
      Zmat <- Zmat.list[[k]]
      Zsmat <- Zsmat.list[[k]]
      phihat <- Estphi_MAR(Zmat,deltavec)
      pivec <- 1/(1+exp( -c(Zsmat%*%phihat) ))
      pimat <- cbind(pimat,pivec)
    }else{
      Zmat <- Zmat.list[[k]]
      Zsmat <- Zsmat.list[[k]]
      phihat <- Estphi(Bmat,Zsmat,deltavec)
      q <- length(phihat)
      est_adj <- function(adj){
        res_int <- integrate(function(uvec){
          sapply(uvec,function(u){
            sum(1/(1+exp(-adj-c(Zmat%*%phihat[-q])-
                           phihat[q]*( yhat+sqrt(s2)*u ))))*
              dnorm(x=u,mean=0,sd=1)
          })
        },-Inf,Inf)
        return(res_int$value-n1)
      }
      M <- 1
      while(TRUE){
        if( est_adj(M)*est_adj(-M)<0 ){
          break
        }else{
          M <- M+1
        }
      }
      res_adj <- uniroot(est_adj,interval=c(-M,M))
      adj <- res_adj$root
      pivec <- 1/(1+exp( -adj-c(Zsmat%*%phihat) ))
      pimat <- cbind(pimat,pivec)
    }
  }
  ##########################################
  pvec <- EstEL(pimat,Aux,n-n1,1/(10*n))
  betahat <- c(solve(t(Xsmat)%*%diag(pvec)%*%Xsmat)%*%t(Xsmat)%*%diag(pvec)%*%ysvec)
  muhat <- sum(pvec*ysvec)/sum(pvec)
  res_quan <- uniroot(Estquan,interval=range(ysvec),ysvec=ysvec,pvec=pvec,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
MCELmin <- function(nuis,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,Bmat,alpha){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  beta <- nuis
  s2 <- s20
  yhat <- c(Xmat%*%beta)
  Xsmat <- Xmat[deltavec,]
  ysvec <- yvec[deltavec]
  ##########################################
  K <- length(is.mar)
  pimat <- NULL
  for(k in 1:K){
    if(is.mar[k]){
      Zmat <- Zmat.list[[k]]
      Zsmat <- Zsmat.list[[k]]
      phihat <- Estphi_MAR(Zmat,deltavec)
      pivec <- 1/(1+exp( -c(Zsmat%*%phihat) ))
      pimat <- cbind(pimat,pivec)
    }else{
      Zmat <- Zmat.list[[k]]
      Zsmat <- Zsmat.list[[k]]
      phihat <- Estphi(Bmat,Zsmat,deltavec)
      q <- length(phihat)
      est_adj <- function(adj){
        res_int <- integrate(function(uvec){
          sapply(uvec,function(u){
            sum(1/(1+exp(-adj-c(Zmat%*%phihat[-q])-
                           phihat[q]*( yhat+sqrt(s2)*u ))))*
              dnorm(x=u,mean=0,sd=1)
          })
        },-Inf,Inf)
        return(res_int$value-n1)
      }
      M <- 1
      while(TRUE){
        if( est_adj(M)*est_adj(-M)<0 ){
          break
        }else{
          M <- M+1
        }
      }
      res_adj <- uniroot(est_adj,interval=c(-M,M))
      adj <- res_adj$root
      pivec <- 1/(1+exp( -adj-c(Zsmat%*%phihat) ))
      pimat <- cbind(pimat,pivec)
    }
  }
  ##########################################
  pvec <- EstEL(pimat,Aux,n-n1,1/(10*n))
  betahat <- c(solve(t(Xsmat)%*%diag(pvec)%*%Xsmat)%*%t(Xsmat)%*%diag(pvec)%*%ysvec)
  muhat <- sum(pvec*ysvec)/sum(pvec)
  res_quan <- uniroot(Estquan,interval=range(ysvec),ysvec=ysvec,pvec=pvec,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,muhat,quanhat)
  nuishat <- c(betahat)
  return(list(thetahat=thetahat,diff=sum((nuis-nuishat)^2)))
}
dif_nuis <- function(nuis,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,Bmat,alpha){
  res_tmp <- MCELmin(nuis,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,Bmat,alpha)
  return(res_tmp$diff)
}