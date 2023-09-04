rm(list=ls())
library(dplyr)
library(ggplot2)
library(reshape)
library(gridExtra)
wd <- "SET_WD"
wd <- "~/Desktop/reproducibility_materials/"
##############################################################################
n_vec <- c(500,2000)
method_vec <- c("CCA","CEL1","CEL2","CEL3","MCELorc","MCELmin")
pn_vec <- c("b0","b1","b2","b3","mu")
nmethod <- length(method_vec)
nn <- length(n_vec)
PS_vec <- 1:3

RES_total <- NULL
for(n in n_vec){
  for(PS in PS_vec){
    for(method in method_vec){
      load(paste0(wd,"output/res_sim2/res2_",method,PS,n/500,".RData"))
      RES$theta <- as.numeric(RES$theta)
      RES$thetahat <- as.numeric(RES$thetahat)
      RES_total <- rbind(RES_total,RES)
    }
  }
}
##############################################################################
# Table S1
tt <- NULL
for(PSt in PS_vec){
  tt_PS <- NULL
  for(nt in n_vec){
    tt_n <- cbind(c(paste0("\\multirow{",nmethod,"}{*}{",nt,"}"),rep("",nmethod-1)),
                  method_vec)
    for(pn in pn_vec){
      RES_tmp <- RES_total %>%
        filter(PS==PSt,n==nt,parname==pn) %>%
        mutate(method=factor(method,levels=method_vec)) %>%
        mutate(bias=thetahat-theta) %>%
        group_by(method) %>%
        summarise(Bias=mean(bias),SD=sd(bias)) %>%
        mutate(Bias=sprintf("%7.4f",Bias)) %>%
        mutate(SD=sprintf("%.4f",SD)) %>%
        as.matrix()
      tt_n <- cbind(tt_n,rep("",nmethod),RES_tmp[,-1])
    }
    tt_PS <- rbind(tt_PS,tt_n)
  }
  tt_PS <- cbind(c(paste0("\\multirow{",nmethod*nn,"}{*}{RM",PSt,"}"),rep("",nmethod*nn-1)),
                 tt_PS)
  tt <- rbind(tt,tt_PS)
}
tt <- apply(tt,1,function(v){paste(v,collapse=" & ")})
tt <- paste(tt,collapse=" \\\\\n")
cat(tt,file=paste0(wd,"output/summary/sim2_table.txt"))
##############################################################################
# Figure S1
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

title_vec <- paste0("(",letters[1:(length(n_vec)*length(PS_vec))],") ",
                    "RM",rep(PS_vec,each=length(n_vec)),
                    ", n=",rep(n_vec,length(PS_vec)))
title_vec <- matrix(title_vec,nrow=length(PS_vec),ncol=length(n_vec),byrow=T)
bp_total <- list()
for(i in 1:length(PS_vec)){
  for(j in 1:length(n_vec)){
    PSt <- PS_vec[i]
    nt <- n_vec[j]
    RES_tmp <- RES_total %>%
      filter(PS==PSt,parname%in%pn_vec,n==nt) %>%
      mutate(method = factor(method,levels=method_vec)) %>%
      mutate(bias = thetahat-theta) %>%
      mutate(PS=paste0("PS",PS)) %>%
      mutate(parname = factor(parname,levels=pn_vec))
    levels(RES_tmp$parname) <- c("beta0","beta1","beta2","beta3","mean")
    bp <- ggplot(data=RES_tmp) +
      geom_boxplot(aes(x=parname,y=bias,fill=method)) +
      geom_abline(slope = 0,intercept = 0,linetype="dashed") +
      theme_bw() + 
      labs(caption=title_vec[i,j]) + 
      theme(plot.caption = element_text(hjust=0.5,size=15)) +
      theme(legend.text = element_text(size=15),
            legend.position="bottom",
            legend.direction="horizontal") +
      guides(fill = guide_legend(nrow = 1)) +
      xlab("") + ylab("Estimates minus truth")
    bp_legend <- get_legend(bp)
    bp <- bp + theme(legend.position="none")
    bp_total <- append(bp_total,list(bp))
  }
}
bp_total <- append(bp_total,list(bp_legend))
lay_mat <- rbind(matrix(1:(length(n_vec)*length(PS_vec)),length(PS_vec),length(n_vec),byrow=T),
                 length(n_vec)*length(PS_vec)+1)
ggsave(filename=paste0(wd,"output/summary/sim2_boxplot.pdf"),
       plot=arrangeGrob(grobs=bp_total,nrow=length(PS_vec)+1,
                        ncol=length(n_vec),layout_matrix=lay_mat,
                        widths=c(rep(5,length(n_vec))),
                        heights=c(rep(3,length(PS_vec)),1)),
       width=5*length(n_vec),height=3*length(PS_vec)+1,device="pdf")
