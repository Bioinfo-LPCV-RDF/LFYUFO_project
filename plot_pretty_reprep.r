## R script to plot nice rep/rep plots for publication
library(ggplot2)
library(cowplot)
# Rscript <WD_PATH> <comparison> <xname> <yname>
# args <- commandArgs(TRUE)
# # library(dplyr, lib.loc="/home/312.3-StrucDev/312.3.1-Commun/R-3.5.0/")
# # library(dplyr)
# print(args)
# WD=args[1]
# wd=args[1]
# label=args[2] # this is the mode: inPeaks or inLibs
# print(label)
# 
# print(WD)


wd<-"//oslo/share_nfs/312-Projets_PCV/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results/comparisons/comp_reprep_2401/" #LFYUFO_cf_reprep
setwd(wd)

# comparison<-"LFYUFO_cfa_LFYUFO_cfb"
# 
# yname="LFYUFO_rep1"
# xname="LFYUFO_rep2"

reprep_plot<-function(comparison,xname,yname) {
  df<-read.table(dir(comparison, full.names=T, pattern="RiL_RiP.tsv$"),header=TRUE,sep="\t") ## read tsv table
  print(dim(df))
  limitsplot=c(min(df[,8],df[,7]),max(df[,8],df[,7])) # cols 8-7 correspond to RIP normalization

  p<-ggplot(df,aes(x=as.double(df[,8]),y=as.double(df[,7])))+
    scale_y_log10(limits=limitsplot)+
    scale_x_log10(limits=limitsplot)+
    geom_point(alpha=0.3,col="grey45")+
    geom_abline(intercept = 0,linetype="longdash",size=0.1)+
    xlab(xname)+
    ylab(yname)+
    theme_classic()+coord_fixed()+
    theme(axis.title=element_text(size=9),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          plot.margin = unit(c(1,15,1,1),"pt"))
  ggsave(filename = paste0(comparison,"_RIP_rep_plot.svg"),device="svg",plot=p,
       width=1000,height=1000,unit="px")

  
  ## RIL plot
  limitsplot_RIL=c(min(df[,10],df[,9]),max(df[,10],df[,9])) # cols 10-9 correspond to RIL normalization
  
  p2<-ggplot(df,aes(x=as.double(df[,10]),y=as.double(df[,9])))+
    scale_y_log10(limits=limitsplot_RIL)+
    scale_x_log10(limits=limitsplot_RIL)+
    geom_point(alpha=0.3,col="grey45")+
    geom_abline(intercept = 0,linetype="longdash",size=0.1)+
    xlab(xname)+
    ylab(yname)+
    theme_classic()+coord_fixed()+
    theme(axis.title=element_text(size=9),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          plot.margin = unit(c(1,15,1,1),"pt"))
  ggsave(filename = paste0(comparison,"_RIL_rep_plot.svg"),device="svg",plot=p2,
         width=1000,height=1000,unit="px")
  
  return(p2)
}


## LFYUFO_cf:
rep1_rep2<-reprep_plot("LFYUFO_cfa_LFYUFO_cfb",yname=expression(paste("LFY+UFO rep1 (",Log[10]," RPKM)")),
                       xname=expression(paste("LFY+UFO rep2 (",Log[10]," RPKM)")))

rep1_rep3<-reprep_plot("LFYUFO_cfa_LFYUFO_cfc",yname=expression(paste("LFY+UFO rep1 (",Log[10]," RPKM)")),
                       xname=expression(paste("LFY+UFO rep3 (",Log[10]," RPKM)")))

rep2_rep3<-reprep_plot("LFYUFO_cfb_LFYUFO_cfc",yname=expression(paste("LFY+UFO rep2 (",Log[10]," RPKM)")),
                       xname=expression(paste("LFY+UFO rep3 (",Log[10]," RPKM)")))


### plot fig S3
a<-ggplot() + theme_void()


rep_plots<-plot_grid(rep1_rep2,rep1_rep3,rep2_rep3,ncol = 3,rel_widths = 1,rel_heights = 1)
ggsave("rep_plots_cowplot.png",rep_plots,device="png",width=3200, height=1200,units="px",dpi=300)
ggsave("rep_plots_cowplot.svg",rep_plots,device="svg",width=3200, height=1200,units="px",dpi=300)
# ggsave("figure_2_withWB_cowplot.svg",f2_cowplot,device="svg",width=3200, height=2200,units="px")

# fig_S3<-ggdraw() +
#   draw_plot(a) +
#   
#   draw_plot(rep1_rep2) +
#   draw_plot(rep1_rep3) +
#   draw_plot(rep2_rep3) +
#   
#   draw_image(mLUBS,0.02,0.21,0.3, 0.3) +
#   draw_image(dLUBS,0.02,0.018,0.4, 0.4) +
#   draw_image(LFY,0.161,-0.08,0.304, 0.304) + 
#   draw_plot(ROC_final,0.41,0,0.36, 0.36) + 
#   draw_plot(a) +
#   geom_segment(aes(x=0.182,y=0.14, xend = 0.182, yend = 0.445),size=0.65,color="grey81",linetype=2)+
#   draw_text(matrices, x = c(0.35,0.35,0.1), y = c(0.365,0.25,0.07), hjust = 0,size = Fontsize+4)+
#   theme(plot.background = element_rect(fill="white", color = NA))+
#   draw_plot_label(c("A","B","C","D","E"), c(0.01, 0.45,0.01,0.46,0.8), c(0.995, 0.995, 0.45,0.40,0.45), size = Fontsize+5) #0.495


# tf<-read.table(dir(comparison, full.names=T, pattern=".tsv$"),header=TRUE) ## read tsv table
# head(tf)

