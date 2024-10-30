########################################################
## Functions
########################################################
CCA <- function(Xmat,yvec,deltavec,alpha){
  Xtmp <- cbind(1,Xmat[deltavec,])
  ysvec <- yvec[deltavec]
  betahat <- c(solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%ysvec)
  yshat <- c(Xtmp%*%betahat)
  shat <- sqrt(mean((ysvec-yshat)^2))
  muhat <- mean(ysvec)
  quanhat <- quantile(ysvec,prob=alpha,names=F)
  thetahat <- c(betahat,shat,muhat,quanhat)
  return(thetahat)
}
CEL <- function(Xmat,yvec,deltavec,Aux,is.mar,Zmat,Zsmat,Bmat,alpha){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  Xtmp <- cbind(1,Xmat[deltavec,])
  ysvec <- yvec[deltavec]
  ##########################################
  if(is.mar){
    phihat <- Estphi_MAR(Zmat,deltavec)
  }else{
    phihat <- Estphi(Bmat,Zsmat,deltavec)
  }
  pisvec <- 1/(1+exp(-c(Zsmat%*%phihat)))
  ##########################################
  pimat <- matrix(pisvec,ncol=1)
  res_el <- EstEL0(pimat,Aux,n-n1)
  pvec <- res_el$pvec
  betahat <- c(solve(t(Xtmp)%*%diag(pvec)%*%Xtmp)%*%t(Xtmp)%*%diag(pvec)%*%ysvec)
  yshat <- c(Xtmp%*%betahat)
  shat <- sqrt(sum(pvec*(ysvec-yshat)^2)/sum(pvec))
  muhat <- sum(pvec*ysvec)/sum(pvec)
  res_quan <- uniroot(Estquan,interval=range(ysvec),ysvec=ysvec,pvec=pvec,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,shat,muhat,quanhat)
  return(thetahat)
}
MCEL_int <- function(nuis_in,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,phihat.list){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  p <- length(nuis_in)
  beta <- nuis_in[-p]
  s2 <- nuis_in[p]^2
  yhat <- c(cbind(1,Xmat)%*%beta)
  K <- length(is.mar)
  ##########################################
  pimat <- NULL
  for(k in 1:K){
    if(is.mar[k]){
      Zmat <- Zmat.list[[k]]
      Zsmat <- Zsmat.list[[k]]
      phihat <- phihat.list[[k]]
      pisvec <- 1/(1+exp( -c(Zsmat%*%phihat) ))
      pivec <- 1/(1+exp( -c(Zmat%*%phihat) ))
      pimat <- cbind(pimat,pisvec - mean(pivec))
    }else{
      H <- 1
      Zmat <- Zmat.list[[k]]
      Zsmat <- Zsmat.list[[k]]
      phihat <- phihat.list[[k]]
      q <- length(phihat)
      pisvec <- 1/(1+exp( -c(Zsmat%*%phihat) ))
      eta_tmp <- c(cbind(Zmat,yhat)%*%phihat)
      res_int <- integrate(function(uvec){
        sapply(uvec,function(u){
          pivec <- 1/(1+exp(-eta_tmp-phihat[q]*sqrt(s2)*u ) )
          return( mean( (1-pivec)^H*pivec )*dnorm(x=u,mean=0,sd=1) )
        })
      },-Inf,Inf)
      pimat <- cbind(pimat,pisvec - res_int$value - sum(1-(1-pisvec)^H)/n)
    }
  }
  res_el <- EstEL1(pimat,Aux,n-n1)
  ##########################################
  pvec <- res_el$pvec
  return(list(pvec=pvec))
}
MCEL <- function(nuis_init,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,phihat.list,alpha){
  n1 <- sum(deltavec)
  Xtmp <- cbind(1,Xmat[deltavec,])
  ysvec <- yvec[deltavec]
  nuis_prev <- nuis_init
  res_tmp <- MCEL_int(nuis_prev,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,phihat.list)
  pvec_prev <- res_tmp$pvec
  is_conv <- T
  count <- 0
  while(TRUE){
    count <- count+1
    beta_new <- c(solve(t(Xtmp)%*%diag(pvec_prev)%*%Xtmp)%*%t(Xtmp)%*%diag(pvec_prev)%*%ysvec)
    yshat <- c(Xtmp%*%beta_new)
    s2_new <- sum(pvec_prev*(ysvec-yshat)^2)
    nuis_new <- c(beta_new,sqrt(s2_new))
    res_tmp <- MCEL_int(nuis_new,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,phihat.list)
    pvec_new <- res_tmp$pvec
    diff <- sqrt(sum((nuis_new-nuis_prev)^2))
    if(diff < 1e-6 ){
      break
    }else{
      nuis_prev <- nuis_new
      pvec_prev <- pvec_new
    }
    if(count >= 100){
      is_conv <- F
      break
    }
  }
  ##########################################
  ysvec <- yvec[deltavec]
  pvec <- res_tmp$pvec
  betahat <- c(solve(t(Xtmp)%*%diag(pvec)%*%Xtmp)%*%t(Xtmp)%*%diag(pvec)%*%ysvec)
  yshat <- c(Xtmp%*%betahat)
  shat <- sqrt(sum(pvec*(ysvec-yshat)^2)/sum(pvec))
  muhat <- sum(pvec*ysvec)/sum(pvec)
  res_quan <- uniroot(Estquan,interval=range(ysvec),ysvec=ysvec,pvec=pvec,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,shat,muhat,quanhat)
  return(thetahat)
}