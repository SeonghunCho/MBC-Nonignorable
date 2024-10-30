rm(list=ls())
########################################################
library(parallel)
library(MASS)
library(stats)
library(dplyr)
###############################################################
# Set the number of cores used in parallel computation for bootstrap.
nCore <- 1
###############################################################
wd <- "SET_WD"
source(paste0(wd,"code/functions_common.R"))
source(paste0(wd,"code/functions_case.R"))
###############################################################
data_NHANES <- read.csv(paste0(wd,"data/NHANES.csv"))
method_vec <- c("CCA","CEL1","CEL2","CEL3","MCEL")
alpha_quan <- 0.5
for(gen in 1:2){
  ###############################################################
  df_tmp <- data_NHANES[data_NHANES$GEN == gen,]
  yvec <- as.vector(scale(df_tmp$DXA))
  X1vec <- as.vector(scale(df_tmp$AGE))
  X2vec <- as.vector(scale(df_tmp$BMI))
  deltavec <- df_tmp$DEL
  Xmat <- cbind(X1vec,X2vec)
  Bmat <- cbind(1,X1vec,X2vec)
  n <- nrow(df_tmp)
  n1 <- sum(deltavec)
  ###############################################################
  indBS <- lapply(1:200,function(sim_num){
    set.seed(1234+sim_num)
    ind <- sample(1:n,n,replace=T)
    return(ind)
  })
  indBS <- do.call("cbind",indBS)
  ###############################################################
  var_name_vec <- c(paste0("b",0:2),"s","mu","q50")
  var_name_vec <- c(var_name_vec,paste0("var_",var_name_vec))
  ###############################################################
  for(method in method_vec){
    cat("gender = ",gen,", ",method," : ",sep="")
    tic <- Sys.time()
    if(method=="CCA"){
      ##################################
      thetahat <- CCA(Xmat=Xmat,yvec=yvec,deltavec=deltavec,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        return(CCA(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],alpha=alpha_quan))
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }else if(method=="CEL1"){
      ##########################################
      Zmat <- cbind(1,X1vec,X2vec)
      Zsmat <- Zmat[deltavec,]
      Aux <- scale(Xmat,center = T,scale = F)[deltavec,]
      thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,Bmat=Bmat,
                      is.mar=T,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        Zmat <- cbind(1,X1vec,X2vec)[indb,]
        Zsmat <- Zmat[deltavec[indb],]
        Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
        return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                   Bmat=Bmat[indb,],is.mar=T,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }else if(method=="CEL2"){
      ##########################################
      Zmat <- cbind(1,X1vec)
      Zsmat <- cbind(Zmat,yvec)[deltavec,]
      Aux <- scale(Xmat,center = T,scale = F)[deltavec,]
      thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,Bmat=Bmat,
                      is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        Zmat <- cbind(1,X1vec)[indb,]
        Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
        Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
        return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                   Bmat=Bmat[indb,],is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }else if(method=="CEL3"){
      ##########################################
      Zmat <- cbind(1,X2vec)
      Zsmat <- cbind(Zmat,yvec)[deltavec,]
      Aux <- scale(Xmat,center = T,scale = F)[deltavec,]
      thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,Bmat=Bmat,
                      is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        Zmat <- cbind(1,X2vec)[indb,]
        Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
        Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
        return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                   Bmat=Bmat[indb,],is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }else if(method=="MCEL"){
      ##########################################
      Zmat.list <- list()
      Zmat.list[[1]] <- cbind(1,X1vec,X2vec)
      Zmat.list[[2]] <- cbind(1,X1vec)
      Zmat.list[[3]] <- cbind(1,X2vec)
      Bmat <- cbind(1,X1vec,X2vec,X1vec^2,X2vec^2)
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
      Xtmp <- cbind(1,Xmat[deltavec,])
      ysvec <- yvec[deltavec]
      beta_init <- as.vector(solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%ysvec)
      yshat <- c(Xtmp%*%beta_init)
      s2_init <- mean((ysvec-yshat)^2)
      nuis_init <- c(beta_init,sqrt(s2_init))
      thetahat <- MCEL(nuis_init=nuis_init,Xmat=Xmat,yvec=yvec,
                       deltavec=deltavec,Aux=Aux,is.mar=is.mar,
                       Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,phihat.list=phihat.list,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        Zmat.list <- list()
        Zmat.list[[1]] <- cbind(1,X1vec,X2vec)[indb,]
        Zmat.list[[2]] <- cbind(1,X1vec)[indb,]
        Zmat.list[[3]] <- cbind(1,X2vec)[indb,]
        Bmat <- cbind(1,X1vec,X2vec,X1vec^2,X2vec^2)[indb,]
        K <- length(Zmat.list)
        Zsmat.list <- list()
        phihat.list <- list()
        for(k in 1:K){
          if(is.mar[k]){
            Zsmat.list[[k]] <- Zmat.list[[k]][deltavec[indb],]
            phihat.list[[k]] <- Estphi_MAR(Zmat.list[[k]],deltavec[indb])
          }else{
            Zsmat.list[[k]] <- cbind(Zmat.list[[k]],yvec[indb])[deltavec[indb],]
            phihat.list[[k]] <- Estphi(Bmat,Zsmat.list[[k]],deltavec[indb])
          }
        }
        Aux <- scale(Xmat[indb,],center = T,scale = F)[deltavec[indb],]
        Xtmpb <- cbind(1,Xmat[indb,])[deltavec[indb],]
        ysvecb <- yvec[indb]
        ysvecb <- ysvecb[deltavec[indb]]
        beta_init <- as.vector(solve(t(Xtmpb)%*%Xtmpb)%*%t(Xtmpb)%*%ysvecb)
        yshat <- c(Xtmpb%*%beta_init)
        s2_init <- mean((ysvecb-yshat)^2)
        nuis_init <- c(beta_init,sqrt(s2_init))
        thetahatBS <- MCEL(nuis_init=nuis_init,Xmat=Xmat[indb,],yvec=yvec[indb],
                           deltavec=deltavec[indb],Aux=Aux,is.mar=is.mar,
                           Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,phihat.list=phihat.list,alpha=alpha_quan)
        return(thetahatBS)
      },mc.cores = 12,mc.preschedule = F)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }
    RES <- c(thetahat,varBS)
    RES <- as.data.frame(t(RES))
    colnames(RES) <- var_name_vec
    RES$Model <- method
    save(list=c("RES","thetaBS"),file=paste0(wd,"output/res_case/",method,gen,".RData"))
    toc <- Sys.time()
    cat(as.numeric(difftime(toc,tic,units = "mins")),"mins\n",sep="")
  }
}