rm(list=ls())
########################################################
library(parallel)
library(MASS)
library(stats)
########################################################
# Set the number of cores used in parallel computation for repetition.
nCore <- 1
########################################################
wd <- "SET_WD"
source(paste0(wd,"code/functions_sim2.R"))
########################################################
beta0 <- c(0.5,0.5,1)
mu0 <- beta0[1]
s20 <- 1/3
s2x <- 1/3
########################################################
method_vec <- c("CCA","CEL1","CEL2","CEL3","MCEL")
n_vec <- c(1000,4000)
PS_vec <- 1:3
gam_vec <- c(0,0.05)
########################################################
for(i1 in 1:length(method_vec)){
  for(i2 in 1:length(n_vec)){
    for(i3 in 1:length(PS_vec)){
      for(i4 in 1:length(gam_vec)){
        method <- method_vec[i1]
        n      <- n_vec[i2]
        PS     <- PS_vec[i3]
        gam    <- gam_vec[i4]
        cat(method,"-",n,"-",PS,"-",gam,"\n",sep="")
        RES <- mclapply(1:1000,function(sim_num){
          set.seed(1234+sim_num)
          ########################################################
          ## Data generation
          ########################################################
          X1vec <- rnorm(n,0,sqrt(s2x))
          X2vec <- rnorm(n,0,sqrt(s2x))
          evec <- rnorm(n,0,sqrt(s20))
          Xmat <- cbind(X1vec,X2vec)
          yvec <- cbind(1,Xmat)%*%beta0+evec + gam*(X1vec^2-s2x)
          if(PS==1){
            phi <- c(0.879,0.5,-0.5)
            Zmat <- cbind(1,X1vec,X2vec)
            etavec <- c(Zmat%*%phi)
          }else if(PS==2){
            phi <- c(1.004,-0.5,-0.25)
            Zmat <- cbind(1,X1vec,yvec)
            etavec <- c(Zmat%*%phi)
          }else if(PS==3){
            phi <- c(0.761,0.5,0.25)
            Zmat <- cbind(1,X2vec,yvec)
            etavec <- c(Zmat%*%phi)
          }
          Prvec <- 1/(1+exp(-etavec))
          Uvec <- runif(n,0,1)
          deltavec <- (Uvec<=Prvec)
          n1 <- sum(deltavec)
          n0 <- n-n1
          ########################################################
          ## Methods
          ########################################################
          if(method=="CCA"){
            ##########################################
            muhat <- CCA(yvec=yvec,deltavec=deltavec)
          }else if(method=="CEL1"){
            ##########################################
            Zmat <- cbind(1,X1vec,X2vec)
            Zsmat <- Zmat[deltavec,]
            Aux <- scale(Xmat)[deltavec,]
            muhat <- CEL(yvec=yvec,deltavec=deltavec,Aux=Aux,is.mar=T,Zmat=Zmat,Zsmat=Zsmat)
          }else if(method=="CEL2"){
            ##########################################
            Zmat <- cbind(1,X1vec)
            Zsmat <- cbind(Zmat,yvec)[deltavec,]
            Aux <- scale(Xmat)[deltavec,]
            muhat <- CEL(yvec=yvec,deltavec=deltavec,Aux=Aux,is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=cbind(1,X1vec,X2vec))
          }else if(method=="CEL3"){
            ##########################################
            Zmat <- cbind(1,X2vec)
            Zsmat <- cbind(Zmat,yvec)[deltavec,]
            Aux <- scale(Xmat)[deltavec,]
            muhat <- CEL(yvec=yvec,deltavec=deltavec,Aux=Aux,is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=cbind(1,X1vec,X2vec))
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
            nuis_init <- as.vector(solve(t(Xtmp)%*%Xtmp)%*%t(Xtmp)%*%ysvec)
            res_tmp <- MCEL(nuis_init=nuis_init,Xmat=Xmat,yvec=yvec,
                            deltavec=deltavec,Aux=Aux,is.mar=is.mar,
                            Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,phihat.list=phihat.list)
            muhat <- res_tmp$muhat
          }
          res_total <- data.frame(n=n,PS=PS,method=method,gam=gam,
                                  sim_num=sim_num,ratio=mean(deltavec),
                                  mu=mu0,muhat=muhat)
          return(res_total)
        },mc.cores = nCore,mc.preschedule = F)
        RES <- do.call("rbind",RES)
        save(list=c("RES"),file=paste0(wd,"output/res_sim2/res_",i1,i2,i3,i4,".RData"))
      }
    }
  }
}