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
source(paste0(wd,"code/functions_sim.R"))
########################################################
method_vec <- c("CCA","CEL1","CEL2","CEL3","MCELorc","MCELmin")
for(method in method_vec){
  for(n in c(500,2000)){
    for(PS in 1:3){
      cat(method,"-",PS,"-",n,"\n",sep="")
      ########################################################
      ## Basic Setting
      ########################################################
      beta0 <- c(0.5,0.3,0.7,0.5)
      s20 <- 1/3
      mux1 <- 0
      mux2 <- 0
      mux3 <- 0
      s2x <- 1/3
      mu0 <- sum(beta0*c(1,mux1,mux2,mux3))
      alpha_quan <- 0.75
      quan0 <- qnorm(p=alpha_quan,mean=mu0,
                     sd=sqrt( sum(beta0[-1]^2)*s2x+s20 ))
      theta0 <- c(beta0,mu0,quan0)
      ########################################################
      ## Main
      ########################################################
      RES <- mclapply(1:1000,function(sim_num){
        set.seed(1234+sim_num)
        ########################################################
        ## Data generation
        ########################################################
        X1vec <- mux1+rnorm(n,0,sqrt(s2x))
        X2vec <- mux2+rnorm(n,0,sqrt(s2x))
        X3vec <- mux3+rnorm(n,0,sqrt(s2x))
        evec <- rnorm(n,0,sqrt(s20))
        Xmat <- cbind(1,X1vec,X2vec,X3vec)
        yvec <- Xmat%*%beta0+evec
        if(PS==1){
          phi <- c(0.879,0.5,0.5)
          Zmat <- cbind(1,X1vec,X2vec)
          etavec <- c(Zmat%*%phi)
        }else if(PS==2){
          phi <- c(0.751,0.5,0.25)
          Zmat <- cbind(1,X1vec,yvec)
          etavec <- c(Zmat%*%phi)
        }else if(PS==3){
          phi <- c(0.757,0.5,0.25)
          Zmat <- cbind(1,X2vec,yvec)
          etavec <- c(Zmat%*%phi)
        }
        Prvec <- 1/(1+exp(-etavec))
        Uvec <- runif(n,0,1)
        deltavec <- (Uvec<=Prvec)
        X1svec <- X1vec[deltavec]
        X2svec <- X2vec[deltavec]
        X3svec <- X3vec[deltavec]
        ysvec <- yvec[deltavec]
        Prsvec <- Prvec[deltavec]
        n1 <- sum(deltavec)
        n0 <- n-n1
        ########################################################
        ## Methods
        ########################################################
        if(method=="CCA"){
          ##########################################
          thetahat <- CCA(Xmat=Xmat,yvec=yvec,deltavec=deltavec,alpha=alpha_quan)
        }else if(method=="CEL1"){
          ##########################################
          Zmat <- cbind(1,X1vec,X2vec)
          Zsmat <- Zmat[deltavec,]
          Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec),X3vec-mean(X3vec))[deltavec,]
          thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                          is.mar=T,Zmat=Zmat,Zsmat=Zsmat,Bmat=cbind(1,X1vec,X2vec),alpha=alpha_quan)
        }else if(method=="CEL2"){
          ##########################################
          Zmat <- cbind(1,X1vec)
          Zsmat <- cbind(Zmat,yvec)[deltavec,]
          Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec),X3vec-mean(X3vec))[deltavec,]
          thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                          is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=cbind(1,X1vec,X2vec),alpha=alpha_quan)
        }else if(method=="CEL3"){
          ##########################################
          Zmat <- cbind(1,X2vec)
          Zsmat <- cbind(Zmat,yvec)[deltavec,]
          Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec),X3vec-mean(X3vec))[deltavec,]
          thetahat <- CEL(Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                          is.mar=F,Zmat=Zmat,Zsmat=Zsmat,Bmat=cbind(1,X1vec,X2vec),alpha=alpha_quan)
        }else if(method=="MCELorc"){
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
          Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec),X3vec-mean(X3vec))[deltavec,]
          thetahat <- MCEL(nuis=beta0,
                           Xmat=Xmat,yvec=yvec,deltavec=deltavec,Aux=Aux,
                           is.mar=is.mar,Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,
                           Bmat=cbind(1,X1vec,X2vec),alpha=alpha_quan)
        }else if(method=="MCELmin"){
          ##########################################
          Zmat.list <- list()
          Zmat.list[[1]] <- cbind(1,X1vec,X2vec,X3vec)
          Zmat.list[[2]] <- cbind(1,X1vec,X2vec)
          Zmat.list[[3]] <- cbind(1,X1vec,X3vec)
          Zmat.list[[4]] <- cbind(1,X2vec,X3vec)
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
          Aux <- cbind(X1vec-mean(X1vec),X2vec-mean(X2vec),X3vec-mean(X3vec))[deltavec,]
          Xsmat <- cbind(1,X1svec,X2svec,X3svec)
          nuis_init <- as.vector(solve(t(Xsmat)%*%Xsmat)%*%t(Xsmat)%*%ysvec)
          res_opt <- optim(par=nuis_init,fn=dif_nuis,
                           Xmat=Xmat,yvec=yvec,
                           deltavec=deltavec,Aux=Aux,is.mar=is.mar,
                           Bmat=cbind(1,X1vec,X2vec,X3vec),
                           Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan)
          nuishat1 <- res_opt$par
          res_MCELmin <- MCELmin(nuis=nuishat1,Xmat=Xmat,yvec=yvec,
                                 deltavec=deltavec,Aux=Aux,is.mar=is.mar,
                                 Bmat=cbind(1,X1vec,X2vec,X3vec),
                                 Zmat.list=Zmat.list,Zsmat.list=Zsmat.list,alpha=alpha_quan)
          thetahat <- res_MCELmin$thetahat
        }
        res_total <- data.frame(n=n,PS=PS,
                                sim_num=sim_num,
                                method=method,
                                parname=c(paste0("b",0:3),"mu","quan"),
                                theta=theta0,thetahat=thetahat)
        return(res_total)
      },mc.cores = nCore,mc.preschedule = F)
      RES <- do.call("rbind",RES)
      save(list=c("RES"),file=paste0(wd,"output/res_sim2/res2_",method,PS,n/500,".RData"))
    }
  }
}