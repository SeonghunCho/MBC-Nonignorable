rm(list=ls())
library(dplyr)
library(ggplot2)
library(reshape)
library(gridExtra)
wd <- "SET_WD"
##############################################################################
PS_vec <- 1:4
method_vec <- c("VAL","CCA","CEL1","CEL2","CEL3","MCEL")
nmethod <- length(method_vec)
pn_vec <- c("mu","quan")

RES_total <- NULL
for(PS in PS_vec){
  for(method in method_vec){
    load(paste0(wd,"res_sim1/res1_",method,PS,".RData"))
    RES$theta <- as.numeric(RES$theta)
    RES$thetahat <- as.numeric(RES$thetahat)
    RES_total <- rbind(RES_total,RES)
  }
}
##############################################################################
## Table 1
tt <- NULL
for(PSt in PS_vec){
  tt_PS <- cbind(c(paste0("\\multirow{",nmethod,"}{*}{RM",PSt,"}"),rep("",nmethod-1)),
                 method_vec)
  for(pn in pn_vec){
    RES_tmp <- RES_total %>%
      filter(PS==PSt,parname==pn) %>%
      mutate(method=factor(method,levels=method_vec)) %>%
      mutate(bias=thetahat-theta) %>%
      group_by(method) %>%
      summarise(Bias=mean(bias),SD=sd(bias),RMSE=sqrt(mean(bias^2)),
                AL=2*qnorm(0.975)*mean(sqrt(varBS)),CP=mean(inc)*100) %>%
      mutate(Bias=sprintf("%7.4f",Bias)) %>%
      mutate(SD=sprintf("%.4f",SD)) %>%
      mutate(RMSE=sprintf("%.4f",RMSE)) %>%
      mutate(AL=sprintf("%.4f",AL)) %>%
      mutate(CP=sprintf("%4.1f",CP)) %>%
      as.matrix()
    tt_PS <- cbind(tt_PS,rep("",nmethod),RES_tmp[,-1])
  }
  tt <- rbind(tt,tt_PS)
}
tt <- apply(tt,1,function(v){paste(v,collapse=" & ")})
tt <- paste(tt,collapse=" \\\\\n")
cat(tt,file=paste0(wd,"summary/sim1_table.txt"))
##############################################################################
## Figure 1
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

PS_vec <- 1:4
title_vec <- paste0("(",c("a","b","c","d"),") RM",1:4)
bp_total <- list()
for(i in 1:length(PS_vec)){
  PSt <- PS_vec[i]
  RES_tmp <- RES_total %>%
    filter(PS==PSt,parname%in%pn_vec) %>%
    mutate(method = factor(method,levels=method_vec)) %>%
    mutate(bias = thetahat-theta) %>%
    mutate(PS=paste0("PS",PS)) %>%
    mutate(parname = factor(parname,levels=pn_vec))
  levels(RES_tmp$parname) <- c("mean","third quartile")
  bp <- ggplot(data=RES_tmp) +
    geom_boxplot(aes(x=parname,y=bias,fill=method)) +
    geom_abline(slope = 0,intercept = 0,linetype="dashed") +
    theme_bw() + 
    labs(caption=title_vec[i]) + 
    theme(plot.caption = element_text(hjust=0.5,size=15)) +
    theme(legend.position="bottom",
          legend.direction="horizontal") +
    guides(fill = guide_legend(nrow = 1)) +
    xlab("") + ylab("Estimates minus truth")
  bp_legend <- get_legend(bp)
  bp <- bp + theme(legend.position="none")
  bp_total <- append(bp_total,list(bp))
}
bp_total <- append(bp_total,list(bp_legend))
lay_mat <- rbind(matrix(1:4,2,2,byrow = T),5)
ggsave(filename=paste0(wd,"summary/sim1_boxplot.pdf"),
       plot=arrangeGrob(grobs=bp_total,nrow=3,ncol=2,layout_matrix=lay_mat,
                        widths=c(5,5),heights=c(3,3,1)),
       width=10,height=6,device="pdf")
##############################################################################