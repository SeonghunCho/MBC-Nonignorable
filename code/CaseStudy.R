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
source(paste0(wd,"code/functions_case.R"))
###############################################################
data_NHANES <- read.csv(paste0(wd,"data/NHANES.csv"))
method_vec <- c("CCA","CEL1","CEL2","CEL3","MCELmin")
for(gen in 1:2){
  ###############################################################
  df_tmp <- data_NHANES[data_NHANES$GEN == gen,]
  yvec <- df_tmp$DXA
  X1vec <- df_tmp$AGE
  X2vec <- df_tmp$BMI
  deltavec <- df_tmp$DEL
  Xmat <- cbind(1,X1vec,X2vec)
  alpha_quan <- 0.5
  n <- nrow(df_tmp)
  n1 <- sum(deltavec)
  ###############################################################
  # There is a computational issue when gen = 1 and bnum = 133.
  indBS <- lapply(1:360,function(sim_num){
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
      Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec))[deltavec,]
      thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,Bmat=Xmat,
                      is.mar=T,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        Zmat <- cbind(1,X1vec,X2vec)[indb,]
        Zsmat <- Zmat[deltavec[indb],]
        Aux <- cbind(X1vec[indb]-mean(X1vec[indb]),
                     X2vec[indb]-mean(X2vec[indb]))[deltavec[indb],]
        return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                   Bmat=Xmat[indb,],is.mar=T,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }else if(method=="CEL2"){
      ##########################################
      Zmat <- cbind(1,X1vec)
      Zsmat <- cbind(Zmat,yvec)[deltavec,]
      Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec))[deltavec,]
      thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,Bmat=Xmat,
                      is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        Zmat <- cbind(1,X1vec)[indb,]
        Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
        Aux <- cbind(X1vec[indb]-mean(X1vec[indb]),
                     X2vec[indb]-mean(X2vec[indb]))[deltavec[indb],]
        return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                   Bmat=Xmat[indb,],is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }else if(method=="CEL3"){
      ##########################################
      Zmat <- cbind(1,X2vec)
      Zsmat <- cbind(Zmat,yvec)[deltavec,]
      Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec))[deltavec,]
      thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,Bmat=Xmat,
                      is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan)
      thetaBS <- mclapply(1:200,function(bnum){
        indb <- indBS[,bnum]
        Zmat <- cbind(1,X1vec)[indb,]
        Zsmat <- cbind(Zmat,yvec[indb])[deltavec[indb],]
        Aux <- cbind(X1vec[indb]-mean(X1vec[indb]),
                     X2vec[indb]-mean(X2vec[indb]))[deltavec[indb],]
        return(CEL(Xmat=Xmat[indb,],yvec=yvec[indb],deltavec=deltavec[indb],Aux=Aux,
                   Bmat=Xmat[indb,],is.mar=F,Zmat=Zmat,Zsmat=Zsmat,alpha=alpha_quan))
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }else if(method=="MCELmin"){
      ##########################################
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
      Xsmat <- Xmat[deltavec,]
      ysvec <- yvec[deltavec]
      beta_init <- as.vector(solve(t(Xsmat)%*%Xsmat)%*%t(Xsmat)%*%ysvec)
      yshat <- c(Xsmat%*%beta_init)
      s2_init <- mean((ysvec-yshat)^2)
      nuis_init <- c(beta_init,s2_init)
      res_opt <- optim(par=nuis_init,fn=dif_nuis,
                       Xmat=Xmat,yvec=yvec,
                       deltavec=deltavec,Aux=Aux,is.mar=is.mar,Bmat=Xmat,
                       Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan)
      nuishat1 <- res_opt$par
      res_MCELmin <- MCELmin(nuis=nuishat1,Xmat=Xmat,yvec=yvec,
                             deltavec=deltavec,Aux=Aux,is.mar=is.mar,Bmat=Xmat,
                             Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan)
      thetahat <- res_MCELmin$thetahat
      if(gen==1){ # There is a computational issue when gen = 1 and bnum = 133.
        ind_vec_bnum <- (1:360)[-(133:144)]
      }else if(gen==2){
        ind_vec_bnum <- 1:360
      }
      thetaBS <- mclapply(ind_vec_bnum,function(bnum){
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
        res_opt <- optim(par=nuis_init,fn=dif_nuis,
                         Xmat=Xmat[indb,],yvec=yvec[indb],
                         deltavec=deltavec[indb],Aux=Aux,is.mar=is.mar,
                         Bmat=Xmat[indb,],
                         Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan)
        nuishat1 <- res_opt$par
        res_MCELmin <- MCELmin(nuis=nuishat1,
                               Xmat=Xmat[indb,],yvec=yvec[indb],
                               deltavec=deltavec[indb],Aux=Aux,is.mar=is.mar,
                               Bmat=Xmat[indb,],
                               Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan)
        thetahat <- res_MCELmin$thetahat
        return(thetahat)
      },mc.preschedule=F,mc.cores=nCore)
      thetaBS <- do.call("rbind",thetaBS)
      varBS <- diag(var(thetaBS))
    }
    RES <- c(thetahat,varBS)
    RES <- as.data.frame(t(RES))
    colnames(RES) <- var_name_vec
    RES$Model <- method
    save(list="RES",file=paste0(wd,"output/res_case/",method,gen,".RData"))
    toc <- Sys.time()
    cat(as.numeric(difftime(toc,tic,units = "mins")),"mins\n",sep="")
  }
}