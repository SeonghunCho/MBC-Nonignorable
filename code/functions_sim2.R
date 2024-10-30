########################################################
## Functions
########################################################
Estphi_MAR <- function(Zmat,deltavec){
  ######################################################
  mllog <- function( x, eps, M, der=0 ){
    # minus log and its first der derivatives, on  eps < x < M
    # 4th order Taylor approx to left of eps and right of M
    # der = 0 or 1 or 2
    # 4th order is lowest that gives self concordance
    
    if( missing(M) )
      M = Inf
    if( eps>M )
      stop("Thresholds out of order")
    
    lo = x < eps
    hi = x > M
    md = (!lo) & (!hi)
    
    # Coefficients for 4th order Taylor approx below eps
    coefs      = rep(0,5)
    coefs[1]   = -log(eps)
    coefs[2:5] = (-eps)^-(1:4)/(1:4)
    
    # Coefficients for 4th order Taylor approx above M
    Coefs      = rep(0,5)
    Coefs[1]   = -log(M)
    Coefs[2:5] = (-M)^-(1:4)/(1:4)
    
    # degree 4 polynomial approx to log
    h = function(y,cvals){ # cvals are coefs at eps, Coefs at M
      # sum c[t+1] y^t
      tee = 1:4
      ans = y*0
      ans = ans + cvals[1]
      for( j in tee )
        ans = ans + y^j*cvals[j+1]
      ans
    }
    
    # first derivative of h at y, from approx at pt
    hp = function(y,pt){
      tee = 0:3
      ans = y*0
      for( j in tee )
        ans = ans + (-y/pt)^j
      ans = ans * (-pt)^-1
      ans
    }
    
    # second derivative of h at y, from approx at pt
    hpp = function(y,pt){
      tee = 0:2
      ans = y*0
      for( j in tee )
        ans = ans + (j+1) * (-y/pt)^j
      ans = ans *(-pt)^-2
      ans
    }
    
    # function value
    f      = x*0
    f[lo]  = h( x[lo]-eps, coefs )
    f[hi]  = h( x[hi]-M,   Coefs )
    f[md]  = -log(x[md])
    
    if( der<1 )return(cbind(f))
    
    # first derivative
    fp     = x*0
    fp[lo] = hp( x[lo]-eps, eps )
    fp[hi] = hp( x[hi]-M, M )
    fp[md] = -1/x[md]
    
    if( der<2 )return(cbind(f,fp))
    
    # second derivative
    fpp     = x*0
    fpp[lo] = hpp( x[lo]-eps, eps )
    fpp[hi] = hpp( x[hi]-M, M )
    fpp[md] = 1/x[md]^2
    
    return( cbind(f,fp,fpp) )
    # End of mllog()
  }
  svdlm <- function(X,y){
    # Linear model regression coefficient via SVD
    
    # Tolerances for generalized inverse via SVD
    RELTOL = 1e-9
    ABSTOL = 1e-100
    
    # Get Xplus = generalized inverse of X
    # If svd algorithm failures are encountered
    # it sometimes helps to try svd(t(X)) and
    # translate back. First check to ensure that
    # X does not contain NaN or Inf or -Inf.
    svdX     = svd(X)
    d        = svdX$d
    lo       = d < (RELTOL * max(d) + ABSTOL)
    dinv     = 1/d
    dinv[lo] = 0
    Xplus    = svdX$v %*% diag(dinv,nrow=length(dinv)) %*% t(svdX$u)
    # taking care with diag when dinv is 1x1
    # to avoid getting the identity matrix of
    # size floor(dinv)
    
    # Return X^+ y
    Xplus %*% matrix(y,ncol=1)
  }
  ######################################################
  ALPHA <- 0.3
  BETA <- 0.8
  BACKEPS <- 0
  eps <- 1e-8
  M <- Inf
  maxit <- 100
  ######################################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  n0 <- n-n1
  phihat <- c(log(n1/n0),rep(0,ncol(Zmat)-1))
  
  # Initial f, g
  pivec <- 1/(1+exp(-c(Zmat%*%phihat)))
  oldvals1 <- mllog(x=pivec,eps=eps,M=M,der=2)
  oldvals2 <- mllog(x=1-pivec,eps=eps,M=M,der=2)
  
  fold <- sum( (deltavec*oldvals1[,1]+(1-deltavec)*oldvals2[,1]) )
  gold <- t(Zmat) %*% ( (deltavec*oldvals1[,2]-(1-deltavec)*oldvals2[,2])*pivec*(1-pivec) )
  
  converged <- FALSE
  iter      <- 0
  while( !converged ){
    iter <- iter + 1
    
    # Get Newton Step
    zt <- t(Zmat) %*% 
      diag( (deltavec*oldvals1[,3]+(1-deltavec)*oldvals2[,3])*pivec^2*(1-pivec)^2+
              (deltavec*oldvals1[,2]-(1-deltavec)*oldvals2[,2])*(1-2*pivec)*pivec*(1-pivec) ) %*% Zmat
    yt <- gold
    step <- -svdlm(zt,yt)  #  more reliable than step = -lm( yt~zt-1 )$coef
    
    backtrack <- FALSE
    s <- 1   # usually called t, but R uses t for transpose
    while( !backtrack ){
      pivecnew <- 1/(1+exp(-c(Zmat%*%(phihat + s*step))))
      newvals1 <- mllog(x=pivecnew,eps=eps,M=M,der=2)
      newvals2 <- mllog(x=1-pivecnew,eps=eps,M=M,der=2)
      fnew <- sum( (deltavec*newvals1[,1]+(1-deltavec)*newvals2[,1]) )
      targ <- fold + ALPHA * s * sum( gold*step ) + BACKEPS # (BACKEPS for roundoff, should not be needed)
      if(  fnew <= targ ){
        # backtracking has converged
        backtrack <- TRUE
        pivec     <- pivecnew
        oldvals1  <- newvals1
        oldvals2  <- newvals2
        fold      <- fnew
        gold      <- t(Zmat) %*% ( (deltavec*oldvals1[,2]-(1-deltavec)*oldvals2[,2])*pivec*(1-pivec) )
        # take the step
        phihat       <- phihat + s*step
      }else{
        s <- s * BETA
      }
    }
    
    # Newton decrement and gradient norm
    ndec     <- sqrt( sum( (step*gold)^2 ) )
    gradnorm <- sqrt( sum(gold^2))
    
    # print(c(fold,gradnorm,ndec,lam))
    
    converged <- ( ndec^2 <= 1e-8)
    if( iter > maxit )break
  }
  return(c(phihat))
}
Estphi <- function(Bmat,Zsmat,deltavec){
  ######################################################
  mreciprocal <- function( x, eps, M, der=0 ){
    # minus log and its first der derivatives, on  eps < x < M
    # 4th order Taylor approx to left of eps and right of M
    # der = 0 or 1 or 2
    # 4th order is lowest that gives self concordance
    
    if( missing(M) )
      M = Inf
    if( eps>M )
      stop("Thresholds out of order")
    
    lo = x < eps
    hi = x > M
    md = (!lo) & (!hi)
    
    # Coefficients for 4th order Taylor approx below eps
    coefs      = -(-1/eps)^(1:5)
    
    # Coefficients for 4th order Taylor approx above M
    Coefs      = -(-1/M)^(1:5)
    
    # degree 4 polynomial approx to log
    h = function(y,cvals){ # cvals are coefs at eps, Coefs at M
      # sum c[t+1] y^t
      tee = 1:4
      ans = y*0
      ans = ans + cvals[1]
      for( j in tee )
        ans = ans + y^j*cvals[j+1]
      ans
    }
    
    # first derivative of h at y, from approx at pt
    hp = function(y,pt){
      tee = 0:3
      ans = y*0
      for( j in tee )
        ans = ans + (-y/pt)^j
      ans = ans * (-pt)^-1
      ans
    }
    
    # second derivative of h at y, from approx at pt
    hpp = function(y,pt){
      tee = 0:2
      ans = y*0
      for( j in tee )
        ans = ans + (j+1) * (-y/pt)^j
      ans = ans *(-pt)^-2
      ans
    }
    
    # function value
    f      = x*0
    f[lo]  = h( x[lo]-eps, coefs )
    f[hi]  = h( x[hi]-M,   Coefs )
    f[md]  = 1/x[md]
    
    if( der<1 )return(cbind(f))
    
    # first derivative
    fp     = x*0
    fp[lo] = hp( x[lo]-eps, eps )
    fp[hi] = hp( x[hi]-M, M )
    fp[md] = -1/x[md]^2
    
    if( der<2 )return(cbind(f,fp))
    
    # second derivative
    fpp     = x*0
    fpp[lo] = hpp( x[lo]-eps, eps )
    fpp[hi] = hpp( x[hi]-M, M )
    fpp[md] = 2/x[md]^3
    
    return( cbind(f,fp,fpp) )
  }
  svdlm <- function(X,y){
    # Linear model regression coefficient via SVD
    
    # Tolerances for generalized inverse via SVD
    RELTOL = 1e-9
    ABSTOL = 1e-100
    
    # Get Xplus = generalized inverse of X
    # If svd algorithm failures are encountered
    # it sometimes helps to try svd(t(X)) and
    # translate back. First check to ensure that
    # X does not contain NaN or Inf or -Inf.
    svdX     = svd(X)
    d        = svdX$d
    lo       = d < (RELTOL * max(d) + ABSTOL)
    dinv     = 1/d
    dinv[lo] = 0
    Xplus    = svdX$v %*% diag(dinv,nrow=length(dinv)) %*% t(svdX$u)
    # taking care with diag when dinv is 1x1
    # to avoid getting the identity matrix of
    # size floor(dinv)
    
    # Return X^+ y
    Xplus %*% matrix(y,ncol=1)
  }
  ######################################################
  ALPHA <- 0.3
  BETA <- 0.8
  # BETA <- 0.5
  BACKEPS <- 0
  eps <- 1e-8
  M <- Inf
  maxit <- 50
  ######################################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  n0 <- n-n1
  Bsmat <- Bmat[deltavec,]
  
  phihat <- c(log(n1/n0),rep(0,ncol(Zsmat)-1))
  q <- length(phihat)
  
  # Initial f, g
  pisvec <- 1/(1+exp(-c(Zsmat%*%phihat)))
  oldvals <- mreciprocal(x=pisvec,eps=eps,M=M,der=2)
  
  fold <- 0.5*sum((t(Bsmat)%*%oldvals[,1] - t(Bmat)%*%rep(1,n))^2)
  gold <- (t(Zsmat)%*%diag(oldvals[,2]*(1-pisvec)*pisvec)%*%Bsmat) %*% 
    (t(Bsmat)%*%oldvals[,1] - t(Bmat)%*%rep(1,n))
  
  converged <- FALSE
  iter      <- 0
  while( !converged ){
    iter <- iter + 1
    
    # Get Newton Step
    zt <- (t(Zsmat)%*%diag(oldvals[,2]*(1-pisvec)*pisvec)%*%Bsmat) %*%
      t( (t(Zsmat)%*%diag(oldvals[,2]*(1-pisvec)*pisvec)%*%Bsmat) )
    for(j in 1:q){
      for(k in 1:q){
        zt[j,k] <- zt[j,k] + 
          sum(c( t(Bsmat)%*%( oldvals[,3]*(1-pisvec)^2*pisvec^2*Zsmat[,j]*Zsmat[,k]+
                                oldvals[,2]*(1-2*pisvec)*(1-pisvec)*pisvec*Zsmat[,j]*Zsmat[,k] ) ) * 
                c(t(Bsmat)%*%oldvals[,1] - t(Bmat)%*%rep(1,n)) )
      }
    }
    yt <- gold
    step <- -svdlm(zt,yt)  #  more reliable than step = -lm( yt~zt-1 )$coef
    
    backtrack <- FALSE
    s <- 1   # usually called t, but R uses t for transpose
    while( !backtrack ){
      pisvecnew <- 1/(1+exp(-c(Zsmat%*%(phihat + s*step))))
      newvals <- mreciprocal(x=pisvecnew,eps=eps,M=M,der=2)
      fnew <- 0.5*sum((t(Bsmat)%*%newvals[,1] - t(Bmat)%*%rep(1,n))^2)
      targ <- fold + ALPHA * s * sum( gold*step ) + BACKEPS # (BACKEPS for roundoff, should not be needed)
      if(  fnew <= targ ){
        # backtracking has converged
        backtrack <- TRUE
        pisvec    <- pisvecnew
        oldvals   <- newvals
        fold      <- fnew
        gold      <- (t(Zsmat)%*%diag(oldvals[,2]*(1-pisvec)*pisvec)%*%Bsmat) %*% 
          (t(Bsmat)%*%oldvals[,1] - t(Bmat)%*%rep(1,n))
        # take the step
        phihat       <- phihat + s*step
      }else{
        s <- s * BETA
      }
    }
    
    # Newton decrement and gradient norm
    ndec     <- sqrt( sum( (step*gold)^2 ) )
    gradnorm <- sqrt( sum(gold^2))
    
    # print(c(fold,gradnorm,ndec,lam))
    
    converged <- ( ndec^2 <= 1e-8)
    if( iter > maxit )break
  }
  return(c(phihat))
}
########################################################
# using W
EstEL0 <- function(P,A,n0,maxit=100){
  mllog <- function( x, eps, M, der=0 ){
    # minus log and its first der derivatives, on  eps < x < M
    # 4th order Taylor approx to left of eps and right of M
    # der = 0 or 1 or 2
    # 4th order is lowest that gives self concordance
    
    if( missing(M) )
      M = Inf
    if( eps>M )
      stop("Thresholds out of order")
    
    lo = x < eps
    hi = x > M
    md = (!lo) & (!hi)
    
    # Coefficients for 4th order Taylor approx below eps
    coefs      = rep(0,5)
    coefs[1]   = -log(eps)
    coefs[2:5] = (-eps)^-(1:4)/(1:4)
    
    # Coefficients for 4th order Taylor approx above M
    Coefs      = rep(0,5)
    Coefs[1]   = -log(M)
    Coefs[2:5] = (-M)^-(1:4)/(1:4)
    
    # degree 4 polynomial approx to log
    h = function(y,cvals){ # cvals are coefs at eps, Coefs at M
      # sum c[t+1] y^t
      tee = 1:4
      ans = y*0
      ans = ans + cvals[1]
      for( j in tee )
        ans = ans + y^j*cvals[j+1]
      ans
    }
    
    # first derivative of h at y, from approx at pt
    hp = function(y,pt){
      tee = 0:3
      ans = y*0
      for( j in tee )
        ans = ans + (-y/pt)^j
      ans = ans * (-pt)^-1
      ans
    }
    
    # second derivative of h at y, from approx at pt
    hpp = function(y,pt){
      tee = 0:2
      ans = y*0
      for( j in tee )
        ans = ans + (j+1) * (-y/pt)^j
      ans = ans *(-pt)^-2
      ans
    }
    
    # function value
    f      = x*0
    f[lo]  = h( x[lo]-eps, coefs )
    f[hi]  = h( x[hi]-M,   Coefs )
    f[md]  = -log(x[md])
    
    if( der<1 )return(cbind(f))
    
    # first derivative
    fp     = x*0
    fp[lo] = hp( x[lo]-eps, eps )
    fp[hi] = hp( x[hi]-M, M )
    fp[md] = -1/x[md]
    
    if( der<2 )return(cbind(f,fp))
    
    # second derivative
    fpp     = x*0
    fpp[lo] = hpp( x[lo]-eps, eps )
    fpp[hi] = hpp( x[hi]-M, M )
    fpp[md] = 1/x[md]^2
    
    return( cbind(f,fp,fpp) )
    # End of mllog()
  }
  svdlm <- function(X,y){
    # Linear model regression coefficient via SVD
    
    # Tolerances for generalized inverse via SVD
    RELTOL = 1e-9
    ABSTOL = 1e-100
    
    # Get Xplus = generalized inverse of X
    # If svd algorithm failures are encountered
    # it sometimes helps to try svd(t(X)) and
    # translate back. First check to ensure that
    # X does not contain NaN or Inf or -Inf.
    svdX     = svd(X)
    d        = svdX$d
    lo       = d < (RELTOL * max(d) + ABSTOL)
    dinv     = 1/d
    dinv[lo] = 0
    Xplus    = svdX$v %*% diag(dinv,nrow=length(dinv)) %*% t(svdX$u)
    # taking care with diag when dinv is 1x1
    # to avoid getting the identity matrix of
    # size floor(dinv)
    
    # Return X^+ y
    Xplus %*% matrix(y,ncol=1)
  }
  
  ALPHA <- 0.3
  BETA <- 0.8
  BACKEPS <- 0
  
  P <- as.matrix(P)
  n1 <- nrow(P)
  K <- ncol(P)
  if(is.null(A)){
    r <- 0
  }else{
    A <- as.matrix(A)
    r <- ncol(A)
  }
  
  n <- n1+n0
  z <- cbind(P-1.0,A)
  
  # eps <- 1/n
  eps <- 1e-8
  M <- Inf
  
  lam <- c( rep(n/(K*n1),K) , rep(0,r) )
  init1 <- mllog( (n/n1)+z%*%lam, eps=eps, M=M, der=2 )
  init2 <- mllog( sum(lam[1:K]), eps=eps, M=M, der=2 )
  
  # Initial f, g
  fold <- sum(init1[,1]) + n0*init2[,1]
  gold <- apply( z * init1[,2],2,sum )
  gold[1:K] <- gold[1:K] + n0*init2[2]
  
  converged <- FALSE
  iter      <- 0
  oldvals1  <- init1
  oldvals2  <- init2
  while( !converged ){
    iter <- iter + 1
    
    # Get Newton Step
    rootllpp <- sqrt(oldvals1[,3])  # sqrt 2nd deriv of -llog lik
    rootllW <- sqrt(n0*oldvals2[3])
    zt <- z
    for( j in 1:(K+r) )
      zt[,j] <- zt[,j] * rootllpp
    zt <- rbind(zt,c(rep(rootllW,K),rep(0.0,r)))
    yt   <- c(oldvals1[,2] / rootllpp, (n0*oldvals2[2])/rootllW)
    step <- -svdlm(zt,yt)  #  more reliable than step = -lm( yt~zt-1 )$coef
    
    backtrack <- FALSE
    s <- 1   # usually called t, but R uses t for transpose
    while( !backtrack ){
      newvals1 <- mllog( (n/n1)+z%*%(lam+s*step),eps=eps,M=M,der=2 )
      newvals2 <- mllog( sum((lam+s*step)[1:K]),eps=eps,M=M,der=2 )
      fnew <- sum(newvals1[,1]) + n0*newvals2[,1]
      targ <- fold + ALPHA * s * sum( gold*step ) + BACKEPS # (BACKEPS for roundoff, should not be needed)
      if(  fnew <= targ ){
        # backtracking has converged
        backtrack <- TRUE
        oldvals1  <- newvals1
        oldvals2  <- newvals2
        fold      <- fnew
        gold      <- apply( z * oldvals1[,2],2,sum )
        gold[1:K] <- gold[1:K] + n0*oldvals2[,2]
        # take the step
        lam       <- lam + s*step
      }else{
        s <- s * BETA
      }
    }
    
    # Newton decrement and gradient norm
    ndec     <- sqrt( sum( (step*gold)^2 ) )
    gradnorm <- sqrt( sum(gold^2))
    
    # print(c(fold,gradnorm,ndec,lam))
    
    converged <- ( ndec^2 <= 1e-8)
    if( iter > maxit )break
  }
  
  pvec <- as.vector(1/(n/n1+z%*%lam))
  pvec <- pvec/sum(pvec)
  return(list("pvec"=pvec,"lam"=as.vector(lam),"elval"=fold))
}
# without using W
EstEL1 <- function(P,A,n0,maxit=100){
  mllog <- function( x, eps, M, der=0 ){
    # minus log and its first der derivatives, on  eps < x < M
    # 4th order Taylor approx to left of eps and right of M
    # der = 0 or 1 or 2
    # 4th order is lowest that gives self concordance
    
    if( missing(M) )
      M = Inf
    if( eps>M )
      stop("Thresholds out of order")
    
    lo = x < eps
    hi = x > M
    md = (!lo) & (!hi)
    
    # Coefficients for 4th order Taylor approx below eps
    coefs      = rep(0,5)
    coefs[1]   = -log(eps)
    coefs[2:5] = (-eps)^-(1:4)/(1:4)
    
    # Coefficients for 4th order Taylor approx above M
    Coefs      = rep(0,5)
    Coefs[1]   = -log(M)
    Coefs[2:5] = (-M)^-(1:4)/(1:4)
    
    # degree 4 polynomial approx to log
    h = function(y,cvals){ # cvals are coefs at eps, Coefs at M
      # sum c[t+1] y^t
      tee = 1:4
      ans = y*0
      ans = ans + cvals[1]
      for( j in tee )
        ans = ans + y^j*cvals[j+1]
      ans
    }
    
    # first derivative of h at y, from approx at pt
    hp = function(y,pt){
      tee = 0:3
      ans = y*0
      for( j in tee )
        ans = ans + (-y/pt)^j
      ans = ans * (-pt)^-1
      ans
    }
    
    # second derivative of h at y, from approx at pt
    hpp = function(y,pt){
      tee = 0:2
      ans = y*0
      for( j in tee )
        ans = ans + (j+1) * (-y/pt)^j
      ans = ans *(-pt)^-2
      ans
    }
    
    # function value
    f      = x*0
    f[lo]  = h( x[lo]-eps, coefs )
    f[hi]  = h( x[hi]-M,   Coefs )
    f[md]  = -log(x[md])
    
    if( der<1 )return(cbind(f))
    
    # first derivative
    fp     = x*0
    fp[lo] = hp( x[lo]-eps, eps )
    fp[hi] = hp( x[hi]-M, M )
    fp[md] = -1/x[md]
    
    if( der<2 )return(cbind(f,fp))
    
    # second derivative
    fpp     = x*0
    fpp[lo] = hpp( x[lo]-eps, eps )
    fpp[hi] = hpp( x[hi]-M, M )
    fpp[md] = 1/x[md]^2
    
    return( cbind(f,fp,fpp) )
    # End of mllog()
  }
  svdlm <- function(X,y){
    # Linear model regression coefficient via SVD
    
    # Tolerances for generalized inverse via SVD
    RELTOL = 1e-9
    ABSTOL = 1e-100
    
    # Get Xplus = generalized inverse of X
    # If svd algorithm failures are encountered
    # it sometimes helps to try svd(t(X)) and
    # translate back. First check to ensure that
    # X does not contain NaN or Inf or -Inf.
    svdX     = svd(X)
    d        = svdX$d
    lo       = d < (RELTOL * max(d) + ABSTOL)
    dinv     = 1/d
    dinv[lo] = 0
    Xplus    = svdX$v %*% diag(dinv,nrow=length(dinv)) %*% t(svdX$u)
    # taking care with diag when dinv is 1x1
    # to avoid getting the identity matrix of
    # size floor(dinv)
    
    # Return X^+ y
    Xplus %*% matrix(y,ncol=1)
  }
  
  ALPHA <- 0.3
  BETA <- 0.8
  BACKEPS <- 0
  
  P <- as.matrix(P)
  n1 <- nrow(P)
  K <- ncol(P)
  if(is.null(A)){
    r <- 0
  }else{
    A <- as.matrix(A)
    r <- ncol(A)
  }
  
  n <- n1+n0
  X <- cbind(P,A)
  
  # eps <- 1/n
  eps <- 1e-8
  M <- Inf
  
  lam <- c( rep(n/(K*n1),K) , rep(0,r) )
  
  init <- mllog( 1+X%*%lam, eps=eps, M=M, der=2 )
  
  # Initial f, g
  fold <- sum(init[,1])
  gold <- t(X)%*%init[,2]
  
  converged <- FALSE
  iter      <- 0
  oldvals  <- init
  while( !converged ){
    iter <- iter + 1
    
    # Get Newton Step
    rootllpp <- sqrt(oldvals[,3])  # sqrt 2nd deriv of -llog lik
    zt <- X
    for( j in 1:(K+r) )
      zt[,j] <- zt[,j] * rootllpp
    yt   <- oldvals[,2] / rootllpp
    step <- -svdlm(zt,yt)  #  more reliable than step = -lm( yt~zt-1 )$coef
    
    backtrack <- FALSE
    s <- 1   # usually called t, but R uses t for transpose
    while( !backtrack ){
      newvals <- mllog( 1+X%*%(lam+s*step),eps=eps,M=M,der=2 )
      fnew <- sum(newvals[,1])
      targ <- fold + ALPHA * s * sum( gold*step ) + BACKEPS # (BACKEPS for roundoff, should not be needed)
      if(  fnew <= targ ){
        # backtracking has converged
        backtrack <- TRUE
        oldvals   <- newvals
        fold      <- fnew
        gold      <- t(X)%*%oldvals[,2]
        # take the step
        lam       <- lam + s*step
      }else{
        s <- s * BETA
      }
    }
    
    # Newton decrement and gradient norm
    ndec     <- sqrt( sum( (step*gold)^2 ) )
    gradnorm <- sqrt( sum(gold^2))
    
    # print(c(fold,gradnorm,ndec,lam))
    
    converged <- ( ndec^2 <= 1e-8)
    if( iter > maxit )break
  }
  
  pvec <- as.vector(1/(1+X%*%lam))
  pvec <- pvec/sum(pvec)
  return(list("pvec"=pvec,"lam"=as.vector(lam),"elval"=fold))
}
########################################################
CCA <- function(yvec,deltavec){
  ysvec <- yvec[deltavec]
  muhat <- mean(ysvec)
  return(muhat)
}
CEL <- function(yvec,deltavec,Aux,is.mar,Zmat,Zsmat,Bmat=NULL){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
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
  muhat <- sum(pvec*ysvec)/sum(pvec)
  return(muhat)
}
MCEL_int <- function(nuis_in,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,phihat.list){
  ##########################################
  n <- length(deltavec)
  n1 <- sum(deltavec)
  beta <- nuis_in
  s2 <- 1/3
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
MCEL <- function(nuis_init,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,phihat.list){
  ##########################################
  n1 <- sum(deltavec)
  Xtmp <- cbind(1,Xmat[deltavec,])
  ysvec <- yvec[deltavec]
  ##########################################
  nuis_prev <- nuis_init
  res_tmp <- MCEL_int(nuis_prev,Xmat,yvec,deltavec,Aux,is.mar,Zmat.list,Zsmat.list,phihat.list)
  pvec_prev <- res_tmp$pvec
  is_conv <- T
  count <- 0
  while(TRUE){
    count <- count+1
    nuis_new <- c(solve(t(Xtmp)%*%diag(pvec_prev)%*%Xtmp)%*%t(Xtmp)%*%diag(pvec_prev)%*%ysvec)
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
  pvec <- res_tmp$pvec
  muhat <- sum(pvec*ysvec)/sum(pvec)
  return(list(muhat=muhat,is_conv=is_conv,count=count))
}