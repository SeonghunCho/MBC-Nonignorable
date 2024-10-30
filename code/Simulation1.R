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
source(paste0(wd,"code/functions_sim1.R"))
########################################################
method_vec <- c("CCA","VAL","IMP","CEL1","CEL2","CEL3","MCEL")
PS_vec <- 1:4

for(i1 in 1:length(method_vec)){
  for(i2 in 1:length(PS_vec)){
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
    sd0 <- sqrt( sum(beta0[2:3]^2)*s2x+s20 )
    quan0 <- qnorm(p=alpha_quan,mean=mu0,sd=sd0)
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
    method <- method_vec[i1]
    PS     <- PS_vec[i2]
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
      Xmat <- cbind(X1vec,X2vec)
      yvec <- cbind(1,Xmat)%*%beta0+evec
      if(PS==1){
        phi <- c(0.879,0.5,-0.5)
        Zmat <- cbind(1,X1vec,X2vec)
        etavec <- c(Zmat%*%phi)
      }else if(PS==2){
        phi <- c(-0.114,0.5,0.25)
        Zmat <- cbind(1,X1vec,yvec)
        etavec <- c(Zmat%*%phi)
      }else if(PS==3){
        phi <- c(0.865,0.5,-0.25)
        Zmat <- cbind(1,X2vec,yvec)
        etavec <- c(Zmat%*%phi)
      }else if(PS==4){
        cutvalue <- mu0
        phi1 <- c(0.857,0.5,-0.25)
        phi2 <- c(0.865,0.5,-0.25)
        Zmat1 <- cbind(1,X1vec,yvec)
        Zmat2 <- cbind(1,X2vec,yvec)
        etavec1 <- c(Zmat1%*%phi1)
        etavec2 <- c(Zmat2%*%phi2)
        etavec <- rep(0,n)
        etavec <- etavec1
        etavec[yvec>=cutvalue] <- etavec2[yvec>=cutvalue]
      }
      Prvec <- 1/(1+exp(-etavec))
      Uvec <- runif(n,0,1)
      deltavec <- (Uvec<=Prvec)
      n1 <- sum(deltavec)
      n0 <- n-n1
      Bmat <- cbind(1,X1vec,X2vec)
      ########################################################
      ## Validation sample
      ########################################################
      X1val <- mux1+rnorm(nval,0,sqrt(s2x))
      X2val <- mux2+rnorm(nval,0,sqrt(s2x))
      eval <- rnorm(nval,0,sqrt(s20))
      Xval <- cbind(X1val,X2val)
      yval <- cbind(1,Xval)%*%beta0+eval
      ########################################################
      ## Methods
      ########################################################
      if(method=="VAL"){
        ##########################################
        thetahat <- VAL(Xval=Xval,yval=yval,alpha=alpha_quan)
        thetaBS <- lapply(1:nBS,function(bnum){
          indb <- indBSval[,bnum]
          thetaBS_tmp <- VAL(Xval=Xval[indb,],yval=yval[indb],alpha=alpha_quan)
          return(thetaBS_tmp)
        })
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="CCA"){
        ##########################################
        thetahat <- CCA(Xmat=Xmat,yvec=yvec,deltavec=deltavec,alpha=alpha_quan)
        thetaBS <- lapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          thetaBS_tmp <- CCA(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],alpha=alpha_quan)
          return(thetaBS_tmp)
        })
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="IMP"){
        ##########################################
        nuisval <- VALnuis(Xval=Xval,yval=yval)
        thetahat <- IMP(nuis=nuisval,Xmat=Xmat,yvec=yvec,deltavec=deltavec,alpha=alpha_quan)
        thetaBS <- lapply(1:nBS,function(bnum){
          indb <- indBSval[,bnum]
          nuisval <- VALnuis(Xval=Xval[indb,],yval=yval[indb])
          indb <- indBS[,bnum]
          thetaBS_tmp <- IMP(nuis=nuisval,Xmat=Xmat[indb,],yvec=yvec[indb],
                             deltavec=deltavec[indb],alpha=alpha_quan)
          return(thetaBS_tmp)
        })
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="CEL1"){
        ##########################################
        Zmat <- cbind(1,X1vec,X2vec)
        Zsmat <- Zmat[deltavec,]
        Aux <- scale(Xmat,center = T,scale = F)[deltavec,]
        thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                        is.mar=T,Zmat=Zmat,Zsmat=Zsmat,Bmat=Bmat,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          Zmat <- cbind(1,X1vec,X2vec)[indb,]
          Zsmat <- Zmat[deltavec[indb],]
          Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
          thetaBS_tmp <- CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                             is.mar=T,Zmat=Zmat,Zsmat=Zsmat,Bmat=Bmat[indb,],alpha=alpha_quan)
          return(thetaBS_tmp)
        },mc.cores = nCore,mc.preschedule = FALSE)
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="CEL2"){
        ##########################################
        Zmat <- cbind(1,X1vec)
        Zsmat <- cbind(Zmat,yvec)[deltavec,]
        Aux <- scale(Xmat,center = T,scale = F)[deltavec,]
        thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                        is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=Bmat,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          Zmat <- cbind(1,X1vec)[indb,]
          Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
          Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
          thetaBS_tmp <- CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                             is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=Bmat[indb,],alpha=alpha_quan)
          return(thetaBS_tmp)
        },mc.cores = nCore,mc.preschedule = FALSE)
        thetaBS <- do.call("rbind",thetaBS)
      }else if(method=="CEL3"){
        ##########################################
        Zmat <- cbind(1,X2vec)
        Zsmat <- cbind(Zmat,yvec)[deltavec,]
        Aux <- scale(Xmat,center = T,scale = F)[deltavec,]
        thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                        is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=Bmat,alpha=alpha_quan)
        thetaBS <- mclapply(1:nBS,function(bnum){
          indb <- indBS[,bnum]
          Zmat <- cbind(1,X2vec)[indb,]
          Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
          Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
          thetaBS_tmp <- CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                             is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=Bmat[indb,],alpha=alpha_quan)
          return(thetaBS_tmp)
        },mc.cores = nCore,mc.preschedule = FALSE)
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
        phihat.list <- list()
        for(k in 1:K){
          if(is.mar[k]){
            Zsmat.list[[k]] <- Zmat.list[[k]][deltavec,]
            phihat.list[[k]] <- Estphi_MAR(Zmat.list[[k]],deltavec)
          }else{
            Zsmat.list[[k]] <- cbind(Zmat.list[[k]],yvec)[deltavec,]
            phihat.list[[k]] <- Estphi(Bmat,Zsmat.list[[k]],deltavec)
          }
        }
        Aux <- scale(Xmat,center = T,scale = F)[deltavec,]
        thetahat <- MCEL(nuis=nuisval,Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                         Zsmat.list=Zsmat.list,phihat.list=phihat.list,alpha=alpha_quan)
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
          phihat.list <- list()
          for(k in 1:K){
            if(is.mar[k]){
              Zsmat.list[[k]] <- Zmat.list[[k]][deltavec[indb],]
              phihat.list[[k]] <- Estphi_MAR(Zmat.list[[k]],deltavec[indb])
            }else{
              Zsmat.list[[k]] <- cbind(Zmat.list[[k]],yvec[indb])[deltavec[indb],]
              phihat.list[[k]] <- Estphi(Bmat[indb,],Zsmat.list[[k]],deltavec[indb])
            }
          }
          Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
          thetaBS_tmp <- MCEL(nuis=nuisval,Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                              is.mar=is.mar,Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,
                              phihat.list=phihat.list,alpha=alpha_quan)
          return(thetaBS_tmp)
        },mc.cores = nCore,mc.preschedule = FALSE)
        thetaBS <- do.call("rbind",thetaBS)
      }
      varBS <- diag(var(thetaBS))
      lb <- thetahat - qnorm(0.975)*sqrt(varBS)
      ub <- thetahat + qnorm(0.975)*sqrt(varBS)
      inc <- (lb <= theta0)&(theta0 <= ub)
      res_total <- data.frame(n=n,PS=PS,method=method,
                              sim_num=sim_num,ratio=mean(deltavec),
                              parname=c(paste0("b",0:2),"mu","quan"),
                              theta=theta0,thetahat=thetahat,
                              varBS=varBS,inc=inc)
      RES <- rbind(RES,res_total)
    }
    save(list=c("RES"),file=paste0(wd,"output/res_sim1/res_",i1,i2,".RData"))
  }
}