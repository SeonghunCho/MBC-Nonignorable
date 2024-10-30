########################################################
## Functions
########################################################
CCA <- function(Xmat,yvec,deltavec,alpha){
  Xtmp <- cbind(1,Xmat)[deltavec,]
  ysvec <- yvec[deltavec]
  betahat <- c(solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%ysvec)
  muhat <- mean(ysvec)
  quanhat <- quantile(ysvec,prob=alpha,names=F)
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
VAL <- function(Xval,yval,alpha){
  Xtmp <- cbind(1,Xval)
  betahat <- c(solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%yval)
  muhat <- mean(yval)
  quanhat <- quantile(yval,prob=alpha,names=F)
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
VALnuis <- function(Xval,yval){
  Xtmp <- cbind(1,Xval)
  betahat <- c(solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%yval)
  nuishat <- c(betahat)
  return(nuishat)
}
IMP <- function(nuis,Xmat,yvec,deltavec,alpha){
  ##########################################\
  beta <- nuis
  s2 <- 1/3
  yhat <- c(cbind(1,Xmat)%*%beta)
  ysvec <- yvec[deltavec]
  ##########################################
  betahat <- beta
  muhat <- mean(yhat)
  Estquan_tmp <- function(quan,yhat,alpha){
    return(mean(pnorm((quan-yhat)/sqrt(s2))) - alpha)
  }
  res_quan <- uniroot(Estquan_tmp,interval=range(ysvec)+c(-3,3),yhat=yhat,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
CEL <- function(Xmat,yvec,deltavec,Aux,is.mar,Zmat,Zsmat,Bmat,alpha){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  Xtmp <- cbind(1,Xmat)[deltavec,]
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
  muhat <- sum(pvec*ysvec)/sum(pvec)
  res_quan <- uniroot(Estquan,interval=range(ysvec)+c(-3,3),ysvec=ysvec,pvec=pvec,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}
MCEL <- function(nuis,Xmat,yvec,deltavec,Aux,Zsmat.list,phihat.list,alpha){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  beta <- nuis
  s2 <- 1/3
  yhat <- c(cbind(1,Xmat)%*%beta)
  Xtmp <- cbind(1,Xmat)[deltavec,]
  ysvec <- yvec[deltavec]
  ##########################################
  K <- length(Zsmat.list)
  pimat <- NULL
  for(k in 1:K){
    Zsmat <- Zsmat.list[[k]]
    phihat <- phihat.list[[k]]
    pisvec <- 1/(1+exp( -c(Zsmat%*%phihat) ))
    pimat <- cbind(pimat,pisvec)
  }
  res_el <- EstEL0(pimat,Aux,n-n1)
  ##########################################
  pvec <- res_el$pvec
  betahat <- c(solve(t(Xtmp)%*%diag(pvec)%*%Xtmp)%*%t(Xtmp)%*%diag(pvec)%*%ysvec)
  muhat <- sum(pvec*ysvec)/sum(pvec)
  res_quan <- uniroot(Estquan,interval=range(ysvec),ysvec=ysvec,pvec=pvec,alpha=alpha)
  quanhat <- res_quan$root
  thetahat <- c(betahat,muhat,quanhat)
  return(thetahat)
}