rm(list=ls())
library(dplyr)
library(ggplot2)
library(reshape)
library(gridExtra)
wd <- "SET_WD"
##############################################################################
## Figure S2
data_NHANES <- read.csv(paste0(wd,"data/NHANES.csv"))
data_NHANES$GEN <- factor(data_NHANES$GEN,levels=c(2,1))
levels(data_NHANES$GEN) <- c("Female","Male")
hg <- ggplot(data_NHANES,aes(x=DXA,fill=GEN,color=GEN)) +
  geom_histogram(alpha=0.5, position="identity",binwidth=1) +
  theme_bw()
ggsave(filename=paste0(wd,"output/summary/case_hist_dxa.pdf"),
       plot=hg,width=8,height=3,device="pdf")
##############################################################################
method_vec <- c("CCA","CEL1","CEL2","CEL3","MCEL")
RES_total <- NULL
for(gen in 1:2){
  df_tmp <- data_NHANES[data_NHANES$GEN == gen,]
  yvec <- (scale(df_tmp$DXA))
  X1vec <- (scale(df_tmp$AGE))
  X2vec <- (scale(df_tmp$BMI))
  my <- attr(yvec,'scaled:center')
  sy <- attr(yvec,'scaled:scale')
  mx1 <- attr(X1vec,'scaled:center')
  sx1 <- attr(X1vec,'scaled:scale')
  mx2 <- attr(X2vec,'scaled:center')
  sx2 <- attr(X2vec,'scaled:scale')
  for(method in method_vec){
    load(paste0(wd,"output/res_case/",method,gen,".RData"))
    thetahat <- as.numeric(RES[,c(paste0("b",0:2),"s","mu","q50")])
    
    thetaBS[,2] <- thetaBS[,2]*sy/sx1
    thetaBS[,3] <- thetaBS[,3]*sy/sx2
    thetaBS[,1] <- thetaBS[,1]*sy - thetaBS[,2]*mx1 - thetaBS[,3]*mx2 + my
    thetaBS[,4] <- thetaBS[,4] * sy
    thetaBS[,5] <- thetaBS[,5] * sy + my
    thetaBS[,6] <- thetaBS[,6] * sy + my
    varBS <- diag(var(thetaBS))
    
    thetahat[2] <- thetahat[2]*sy/sx1
    thetahat[3] <- thetahat[3]*sy/sx2
    thetahat[1] <- thetahat[1]*sy - thetahat[2]*mx1 - thetahat[3]*mx2 + my
    thetahat[4] <- thetahat[4] * sy
    thetahat[5] <- thetahat[5] * sy + my
    thetahat[6] <- thetahat[6] * sy + my
    
    var_name_vec <- c(paste0("b",0:2),"s","mu","q50")
    var_name_vec <- c(var_name_vec,paste0("var_",var_name_vec))
    
    RES <- c(thetahat,varBS)
    RES <- as.data.frame(t(RES))
    colnames(RES) <- var_name_vec
    RES$Model <- method
    RES$gen <- gen
    RES_total <- rbind(RES_total,RES)
  }
}
##############################################################################
## Table 2
var_name_vec <- c(paste0("b",0:2),"s","mu","q50")
for(gent in 1:2){
  RES_tmp <- RES_total %>% filter(gen==gent)
  
  tt <- RES_tmp$Model
  for(vn in var_name_vec){
    est <- sprintf("%6.3f",RES_tmp[,vn])
    est_sd <- sprintf("%5.3f",sqrt(RES_tmp[,paste0("var_",vn)]))
    tt <- cbind(tt,"",est,est_sd)
  }
  tt <- apply(tt,1,function(v){paste(v,collapse=" & ")})
  tt <- paste(tt,collapse=" \\\\\n")
  cat(tt,file=paste0(wd,"output/summary/case_table",gent,".txt"))
}
##############################################################################
## Figure S3
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
RES_df <- melt(RES_total[,c(var_name_vec,"Model","gen")],id.vars=c("Model","gen"))
tmp_df <- melt(RES_total[,c(paste0("var_",var_name_vec),"Model","gen")],id.vars=c("Model","gen"))
colnames(RES_df) <- c("method","gen","parname","est")
RES_df$var_est <- tmp_df$value
title_vec <- paste0("(",letters[1:length(var_name_vec)],") ",
                    c("Intercept","Age","BMI","SD of Error","Mean","Median"))
plot_total <- list()
for(i in 1:length(var_name_vec)){
  dl <- 0.4
  vn <- var_name_vec[i]
  RES_tmp <- RES_df %>%
    filter(parname == vn) %>%
    mutate(method = factor(method,levels=method_vec)) %>%
    mutate(gen = factor(gen,levels=1:2)) %>%
    mutate(lb = est-qnorm(0.975)*sqrt(var_est)) %>%
    mutate(ub = est+qnorm(0.975)*sqrt(var_est))
  levels(RES_tmp$gen) <- c("Male","Female")
  plot_tmp <- ggplot(data=RES_tmp,aes(x=gen,y=est,group=method,color=method)) + 
    geom_point(position=position_dodge(dl)) +
    geom_errorbar(aes(ymin=lb,ymax=ub),position="dodge",width=dl) +
    labs(caption=title_vec[i]) + 
    xlab("") + ylab("") + theme_bw() +
    theme(plot.caption = element_text(hjust=0.5,size=15)) +
    theme(legend.text = element_text(size=15),
          legend.position="bottom",
          legend.direction="horizontal") +
    guides(fill = guide_legend(nrow = 1))
  legend_tmp <- get_legend(plot_tmp)
  
  plot_tmp <- plot_tmp + theme(legend.position="none")
  plot_total <- append(plot_total,list(plot_tmp))
}
plot_total <- append(plot_total,list(legend_tmp))
lay_mat <- rbind(matrix(1:length(var_name_vec),ncol=2,byrow=T),length(var_name_vec)+1)
ggsave(filename=paste0(wd,"output/summary/case_errorbarplot.pdf"),
       plot=arrangeGrob(grobs=plot_total,nrow=length(var_name_vec)/2+1,
                        ncol=2,layout_matrix=lay_mat,
                        widths=c(5,5),heights=c(3,3,3,1)),
       width=5*2,height=3*3+1,device="pdf")
##############################################################################