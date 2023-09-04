rm(list=ls())
########################################################
library(parallel)
library(MASS)
library(stats)
########################################################
# Set the number of cores used in parallel computation for bootstrap.
nCore <- 1
########################################################
wd <- "SET_WD"
source(paste0(wd,"code/functions_sim.R"))
########################################################
method_vec <- c("CCA","VAL","CEL1","CEL2","CEL3","MCEL")
for(method in method_vec){
  for(PS in 1:4){
    ########################################################
    ## Basic Setting
    ########################################################
    n <- 2000
    nval <- 0.1*n
    
    beta0 <- c(0.5,1,0.5)
    s20 <- 1/3
    mux1 <- 1
    mux2 <- 1
    s2x <- 1/3
    mu0 <- sum(beta0*c(1,mux1,mux2))
    alpha_quan <- 0.75
    quan0 <- qnorm(p=alpha_quan,mean=mu0,sd=sqrt( sum(beta0[2:3]^2)*s2x+s20 ))
    theta0 <- c(beta0,mu0,quan0)
    ########################################################
    ## Bootstrap index
    ########################################################
    nBS <- 200
    indBS <- lapply(1:nBS,function(sim_num){
      set.seed(1234+sim_num)
      ind <- sample(1:n,n,replace=T)
      return(ind)
    })
    indBS <- do.call("cbind",indBS)
    indBSval <- lapply(1:nBS,function(sim_num){
      set.seed(1234+sim_num)
      ind <- sample(1:nval,nval,replace=T)
      return(ind)
    })
    indBSval <- do.call("cbind",indBSval)
    ########################################################
    ## Main
    ########################################################
    RES <- NULL
    for(sim_num in 1:1000){
      cat(method,"-",PS," : ",sim_num,"th iter\n",sep="")
      set.seed(1234+sim_num)
      ########################################################
      ## Data generation
      ########################################################
      X1vec <- mux1+rnorm(n,0,sqrt(s2x))
      X2vec <- mux2+rnorm(n,0,sqrt(s2x))
      evec <- rnorm(n,0,sqrt(s20))
      Xmat <- cbind(1,X1vec,X2vec)
      yvec <- Xmat%*%beta0+evec
      if(PS==1){
        phi <- c(0.879,0.5,-0.5)
        Zmat <- cbind(1,X1vec,X2vec)
        etavec <- c(Zmat%*%phi)
      }else if(PS==2){
        phi <- c(0.857,0.5,-0.25)
        Zmat <- cbind(1,X1vec,yvec)
        etavec <- c(Zmat%*%phi)
      }else if(PS==3){
        phi <- c(0.865,0.5,-0.25)
        Zmat <- cbind(1,X2vec,yvec)
        etavec <- c(Zmat%*%phi)
      }else if(PS==4){
        cutvalue <- mu0
        phi1 <- c(0.879,0.5,-0.5)
        phi2 <- c(0.857,0.5,-0.25)
        Zmat1 <- cbind(1,X1vec,X2vec)
        Zmat2 <- cbind(1,X1vec,yvec)
        etavec1 <- c(Zmat1%*%phi1)
        etavec2 <- c(Zmat2%*%phi2)
        etavec <- rep(0,n)
        etavec[yvec<cutvalue] <- etavec1[yvec<cutvalue]
        etavec[yvec>=cutvalue] <- etavec2[yvec>=cutvalue]
      }
      Prvec <- 1/(1+exp(-etavec))
      Uvec <- runif(n,0,1)
      deltavec <- (Uvec<=Prvec)
      X1svec <- X1vec[deltavec]
      X2svec <- X2vec[deltavec]
      ysvec <- yvec[deltavec]
      Prsvec <- Prvec[deltavec]
      n1 <- sum(deltavec)
      n0 <- n-n1
      ########################################################
      ## Validation sample
      ########################################################
      X1val <- mux1+rnorm(nval,0,sqrt(s2x))
      X2val <- mux2+rnorm(nval,0,sqrt(s2x))
      eval <- rnorm(nval,0,sqrt(s20))
      Xval <- cbind(1,X1val,X2val)
      yval <- Xval%*%beta0+eval
      ########################################################
      ## Methods
      ########################################################
      if(method=="VAL"){
        ##########################################
        thetahat <- VAL(Xval=Xval,yval=yval,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBSval[,bnum]
          return(VAL(Xval=Xval[indb,],yval=yval[indb],alpha=alpha_quan))
        },mc.preschedule=F,mc.cores=nCore)
        thetaBS <- do.call("rbind",thetaBS)
        
      }else if(method=="CCA"){
        ##########################################
        thetahat <- CCA(Xmat=Xmat,yvec=yvec,deltavec=deltavec,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          return(CCA(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],alpha=alpha_quan))
        },mc.preschedule=F,mc.cores=nCore)
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="CEL1"){
        ##########################################
        Zmat <- cbind(1,X1vec,X2vec)
        Zsmat <- Zmat[deltavec,]
        Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec))[deltavec,]
        thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                        is.mar=T,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          Zmat <- cbind(1,X1vec,X2vec)[indb,]
          Zsmat <- Zmat[deltavec[indb],]
          Aux <- cbind(X1vec[indb]-mean(X1vec[indb]),
                       X2vec[indb]-mean(X2vec[indb]))[deltavec[indb],]
          return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                     is.mar=T,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
        },mc.preschedule=F,mc.cores=nCore)
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="CEL2"){
        ##########################################
        Zmat <- cbind(1,X1vec)
        Zsmat <- cbind(Zmat,yvec)[deltavec,]
        Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec))[deltavec,]
        thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                        is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          Zmat <- cbind(1,X1vec)[indb,]
          Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
          Aux <- cbind(X1vec[indb]-mean(X1vec[indb]),
                       X2vec[indb]-mean(X2vec[indb]))[deltavec[indb],]
          return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                     is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
        },mc.preschedule=F,mc.cores=nCore)
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="CEL3"){
        ##########################################
        Zmat <- cbind(1,X2vec)
        Zsmat <- cbind(Zmat,yvec)[deltavec,]
        Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec))[deltavec,]
        thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                        is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          Zmat <- cbind(1,X2vec)[indb,]
          Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
          Aux <- cbind(X1vec[indb]-mean(X1vec[indb]),
                       X2vec[indb]-mean(X2vec[indb]))[deltavec[indb],]
          return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                     is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
        },mc.preschedule=F,mc.cores=nCore)
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="MCEL"){
        ##########################################
        nuisval <- VALnuis(Xval=Xval,yval=yval)
        Zmat.list <- list()
        Zmat.list[[1]] <- cbind(1,X1vec,X2vec)
        Zmat.list[[2]] <- cbind(1,X1vec)
        Zmat.list[[3]] <- cbind(1,X2vec)
        K <- length(Zmat.list)
        is.mar <- rep(F,K)
        is.mar[1] <- T
        Zsmat.list <- list()
        for(k in 1:K){
          if(is.mar[k]){
            Zsmat.list[[k]] <- Zmat.list[[k]][deltavec,]
          }else{
            Zsmat.list[[k]] <- cbind(Zmat.list[[k]],yvec)[deltavec,]
          }
        }
        Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec))[deltavec,]
        thetahat <- MCEL(nuis=nuisval,
                         Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                         is.mar=is.mar,Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBSval[,bnum]
          nuisval <- VALnuis(Xval=Xval[indb,],yval=yval[indb])
          indb <- indBS[,bnum]
          Zmat.list <- list()
          Zmat.list[[1]] <- cbind(1,X1vec,X2vec)[indb,]
          Zmat.list[[2]] <- cbind(1,X1vec)[indb,]
          Zmat.list[[3]] <- cbind(1,X2vec)[indb,]
          K <- length(Zmat.list)
          Zsmat.list <- list()
          for(k in 1:K){
            if(is.mar[k]){
              Zsmat.list[[k]] <- Zmat.list[[k]][deltavec[indb],]
            }else{
              Zsmat.list[[k]] <- cbind(Zmat.list[[k]],yvec[indb])[deltavec[indb],]
            }
          }
          Aux <- cbind(X1vec[indb]-mean(X1vec[indb]),
                       X2vec[indb]-mean(X2vec[indb]))[deltavec[indb],]
          return(MCEL(nuis=nuisval,
                      Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                      is.mar=is.mar,Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan))
        },mc.preschedule=F,mc.cores=nCore)
        thetaBS <- do.call("rbind",thetaBS)
      }
      varBS <- diag(var(thetaBS))
      lb <- thetahat - qnorm(0.975)*sqrt(varBS)
      ub <- thetahat + qnorm(0.975)*sqrt(varBS)
      inc <- (lb <= theta0)&(theta0 <= ub)
      res_total <- data.frame(n=n,PS=PS,
                              sim_num=sim_num,
                              method=method,
                              parname=c(paste0("b",0:2),"mu","quan"),
                              theta=theta0,thetahat=thetahat,
                              varBS=varBS,inc=inc)
      RES <- rbind(RES,res_total)
    }
    save(list=c("RES"),file=paste0(wd,"output/res_sim1/res1_",method,PS,".RData"))
  }
}