rm(list=ls())
library(dplyr)
library(ggplot2)
library(reshape)
library(gridExtra)
wd <- "SET_WD"
##############################################################################
method_vec <- c("CCA","CEL1","CEL2","CEL3","MCEL")
n_vec <- c(1000,4000)
PS_vec <- 1:3
gam_vec <- c(0,0.05)
nn <- length(n_vec)
nmethod <- length(method_vec)
##############################################################################
RES_total <- NULL
for(i1 in 1:length(method_vec)){
  for(i2 in 1:length(n_vec)){
    for(i3 in 1:length(PS_vec)){
      for(i4 in 1:length(gam_vec)){
        load(paste0(wd,"output/res_sim2/res_",i1,i2,i3,i4,".RData"))
        RES_total <- rbind(RES_total,RES)
      }
    }
  }
}
RES_total <- RES_total  %>%
  mutate(bias = muhat-mu)
##############################################################################
# Table S1
tt <- NULL
for(i2 in 1:length(n_vec)){
  tt_n <- cbind(c(paste0("\\multirow{",nmethod,"}{*}{",n_vec[i2],"}"),rep("",nmethod-1)),
                method_vec)
  for(i3 in 1:length(PS_vec)){
    RES_tmp <- RES_total %>%
      filter(PS==PS_vec[i3],n==n_vec[i2],parname==pn) %>%
      mutate(method=factor(method,levels=method_vec)) %>%
      group_by(method) %>%
      summarise(Bias=mean(bias),SD=sd(bias)) %>%
      mutate(Bias=sprintf("%7.4f",Bias)) %>%
      mutate(SD=sprintf("%.4f",SD)) %>%
      as.matrix()
    tt_n <- cbind(tt_n,rep("",nmethod),RES_tmp[,-1])
    # tt_n <- cbind(tt_n,RES_tmp[,-1])
  }
  tt <- rbind(tt,tt_n)
}
tt <- apply(tt,1,function(v){paste(v,collapse=" & ")})
tt <- paste(tt,collapse=" \\\\\n")
cat(tt,file=paste0(wd,"output/summary/sim2_table.txt"))
##############################################################################
# Figure S1
for(i4 in 1:length(gam_vec)){
  title_vec <- paste0("(",letters[1:(length(n_vec)*length(PS_vec))],") ",
                      "RM",rep(PS_vec,each=length(n_vec)),", n=",rep(n_vec,length(PS_vec)))
  title_vec <- matrix(title_vec,nrow=length(PS_vec),ncol=length(n_vec),byrow=T)
  bp_total <- list()
  for(i3 in 1:length(PS_vec)){
    bias_max <- max(abs(RES_total$bias[ (RES_total$PS==PS_vec[i3]) ]))
    for(i2 in 1:length(n_vec)){
      nt <- n_vec[i2]
      PSt <- PS_vec[i3]
      gamt <- gam_vec[i4]
      RES_tmp <- RES_total %>%
        filter(PS==PSt,n==nt,gam==gamt) %>%
        mutate(method = factor(method,levels=method_vec)) %>%
        mutate(PS=paste0("PS",PS))
      bp <- ggplot(data=RES_tmp) +
        geom_boxplot(aes(x=method,y=bias)) +
        geom_abline(slope = 0,intercept = 0,linetype="dashed") +
        theme_bw() + 
        labs(caption=title_vec[i3,i2]) + 
        theme(plot.caption = element_text(hjust=0.5,size=15)) +
        theme(legend.text = element_text(size=15),
              legend.position="bottom",
              legend.direction="horizontal") +
        guides(fill = guide_legend(nrow = 1)) +
        xlab("") + ylab("Estimates minus truth") +
        ylim(-0.3,0.3)
      bp <- bp + theme(legend.position="none")
      bp_total <- append(bp_total,list(bp))
    }
  }
  lay_mat <- rbind(matrix(1:(length(n_vec)*length(PS_vec)),length(PS_vec),length(n_vec),byrow=T))
  ggsave(filename=paste0(wd,"output/summary/sim2_boxplot",i4,".pdf"),
         plot=arrangeGrob(grobs=bp_total,nrow=length(PS_vec),
                          ncol=length(n_vec),layout_matrix=lay_mat,
                          widths=c(rep(5,length(n_vec))),
                          heights=c(rep(3,length(PS_vec)))),
         width=5*length(n_vec),height=3*length(PS_vec),device="pdf")
}
##############################################################################