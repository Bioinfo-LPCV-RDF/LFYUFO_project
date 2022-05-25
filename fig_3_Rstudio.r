### script to plot figure 3
library(reshape2)
library(plotROC)
library(ggplot2)
library(ggpubr)
library(viridis) #viridislite
library(png) #libpng
library(cowplot)
theme_set(theme_cowplot())

LFYUFO_dir<-"//oslo/share_nfs/312-Projets_PCV/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO"
setwd(paste0(LFYUFO_dir,"/results/comparisons"))

name1<-"LFYamp"
name2<-"LFYUFO_cf"

setwd(paste0("./",name1,"_",name2))

table<-read.table(paste0("table_",name1,"_",name2,"_RiL_RiP.tsv"),header=TRUE)
table$CFC<-(table$LFYUFO_cf_RiP/table$LFYamp_RiP)
table<-table[order(-table$CFC),]


table$CFC_RIL<-(table$LFYUFO_cf_RiL/table$LFYamp_RiL)


## Covplot LFY vs LFYUFO
### change shape of 600 peaks with highest CFC (used for motifs)
twentyperc<-round(dim(table)[1]*0.2)
LFYUFO_spe_vec<-rep("LFYUFO_spe",twentyperc)
other_peaks_vec<-rep("LFY_spe",(dim(table)[1]-twentyperc))
CFC_name<-c(LFYUFO_spe_vec,other_peaks_vec)
table<-cbind(table,CFC_name)


motif_peaks<-rep("training",600)
testing_peaks<-rep("testing",(twentyperc-600))
other_peaks<-rep("other",(dim(table)[1]-twentyperc))
peak_sets<-c(motif_peaks,testing_peaks,other_peaks)
table<-cbind(table,peak_sets)


training_peaks<-rep("training",600) #rep(60,600) 
other_peaks<-rep("other",(dim(table)[1]-600)) #rep(2,(dim(table)[1]-600))
training_set<-c(training_peaks,other_peaks)
table<-cbind(table,training_set)



table$training_set<-factor(table$training_set, levels = c("training","other"))

LFYvsLFYUFO_gradient=c("#ff4c91", "#ff4c60","#ff6a4c", "#ff8a4c", "#ffab4c", "#ffbc4c", "#ffcc4c","#ffdb4c") #"#ff4c4f" "#ff9b4c",
Fontsize=8
limitsplot=c(2,20000)


covplot<-ggplot(data=table,aes(x=LFYUFO_cf_RiP, y=LFYamp_RiP))+
  scale_y_log10(limits=limitsplot)+
  scale_x_log10(limits=limitsplot)+
  geom_point(aes(x=LFYUFO_cf_RiP,y=LFYamp_RiP,col=CFC,shape=training_set,size=training_set),
             alpha=0.3)+
  scale_size_manual(values=c(1.8,2.3),
                    guide="none")+
  scale_colour_gradientn(colours=LFYvsLFYUFO_gradient,
                       name=expression(paste("CFC = ", frac("LFY-UFO","LFY"))),
                       breaks=c(1,5,10,15,20))+
  theme_classic()+
  theme(axis.title = element_text(size=Fontsize+1),
        axis.text = element_text(size=Fontsize),
        legend.title = element_text(size=Fontsize),
        legend.text = element_text(size=Fontsize))+
  xlab(expression(paste("LFY-UFO (",Log[10]," RPKM)")))+
  ylab(expression(paste("LFY (",Log[10]," RPKM)")))+
  scale_shape_manual(values=c(17,19),name="Peaks",
                     labels=c("training set","other"))+
  geom_abline(intercept = 0,linetype="longdash",size=0.1)+
  coord_fixed()+
  guides(color = guide_colourbar(order=1),
         shape = guide_legend(order=2))



####### Cowplot definitive plots #######
## Upload images
dLUBS<-readPNG(source = paste0(LFYUFO_dir,"/results/motifs/dLUBS/LFYUFO_cf/LFYUFO_cf_logo.png")) # new
mLUBS<-readPNG(source = paste0(LFYUFO_dir,"/results/motifs/mLUBS/LFYUFO_cf/LFYUFO_cf_logo.png")) # new
LFY<-readPNG(source = paste0(LFYUFO_dir,"/results/motifs/LFYamp/LFYamp/LFYamp_logo.png")) # new



#################### ROC ####################
### plotROC
#### load external scores files
dLUBS_scores<-read.table(paste0(LFYUFO_dir,"/results/ROCS/LFYUFO_cf_top20CFC/scores/tab_dLUBS.tsv"),header=FALSE)
mLUBS_scores<-read.table(paste0(LFYUFO_dir,"/results/ROCS/LFYUFO_cf_top20CFC/scores/tab_mLUBS.tsv"),header=FALSE)
LFY_scores<-read.table(paste0(LFYUFO_dir,"/results/ROCS/LFYUFO_cf_top20CFC/scores/tab_LFY.tsv"),header=FALSE)


## reshape df
melt_df<-function(df){
  names(df)<-c("pos","neg")
  df$id<-row.names(df)
  melted <- melt(df, id="id")
  melted$D <- ifelse(melted$variable == 'pos',1,0)
  return(melted)
}

dLUBS_melted<-melt_df(dLUBS_scores)
mLUBS_melted<-melt_df(mLUBS_scores)
LFY_melted<-melt_df(LFY_scores)



## melt the three datasets together
full_df<-cbind(dLUBS_melted,mLUBS_melted[,3],LFY_melted[,3])
names(full_df)<-c("id","variable","dLUBS","D","mLUBS","LFY")
full_df<-full_df[,c(1,2,4,3,5,6)]

melted_full_df<-melt_roc(full_df,"D",c("dLUBS","mLUBS","LFY"))



ROC_colors<-c("#09495e","#118ab2","#ff4c91")


melted_full_df$name2 <- factor(melted_full_df$name, levels=c("dLUBS", "mLUBS", "LFY"), labels=c("dLUBS", "mLUBS", "LFY"))

ROC<-ggplot(melted_full_df, aes(d = D, m = M, color = name2,linetype=name2)) +
  geom_roc(n.cuts=0)+
  scale_linetype_manual(values=c("solid","solid","solid"))+
  scale_color_manual(values=ROC_colors)+
  theme_classic()+
  xlab("FPR")+ylab("TPR")+
  
  theme(axis.title = element_text(size=Fontsize+3),
        axis.text = element_text(size=Fontsize+2),
        legend.position = "none")+ coord_fixed()


## add AUC values on graph:
AUCs<-calc_auc(ggroc = ROC)

AUC_dLUBS<-round(AUCs$AUC[AUCs$group==1],digits=2)
AUC_mLUBS<-round(AUCs$AUC[AUCs$group==2],digits=2)
AUC_LFY<-format(round(AUCs$AUC[AUCs$group==3],digits=2), nsmall = 2)


AUC_on_plot<-c(paste0("dLUBS: ",AUC_dLUBS),
               paste0("mLUBS: ",AUC_mLUBS),
               paste0("LFY: ",AUC_LFY))


ROC_final<-ggplot(melted_full_df, aes(d = D, m = M, color = name2,linetype=name2)) +
  geom_roc(n.cuts=0)+
  scale_linetype_manual(values=c("solid","solid","solid"))+
  scale_color_manual(values=ROC_colors)+
  theme_classic()+
  xlab("FPR")+ylab("TPR")+
  
  annotate(geom="text", x = 0.9,y=c(0.25,0.18,0.11), label=AUC_on_plot,
  color=ROC_colors,hjust=1,size=Fontsize-4)+
  
  theme(axis.title = element_text(size=Fontsize),
        axis.text = element_text(size=Fontsize),
        legend.position = "none")+ coord_fixed()




## cowplot multiplot figure
a<-ggplot() + theme_void()
matrices<-c("mLUBS","dLUBS","LFY\ncanonical")

f2_cowplot<-ggdraw() +
  draw_plot(a) +
  
  draw_plot(covplot,0.44,0.45, 0.538, 0.539) +
  
  draw_image(mLUBS,0.02,0.21,0.29, 0.29) +
  draw_image(dLUBS,0.031,0.027,0.38, 0.38) +
  draw_image(LFY,0.156,-0.07,0.3, 0.3) + 
  draw_plot(ROC_final,0.4,0,0.36, 0.36) + 
  draw_plot(a) +
  geom_segment(aes(x=0.178,xend = 0.178,y=0.145,yend = 0.418),size=0.65,color="grey81",linetype=2)+
  draw_text(matrices, x = c(0.35,0.35,0.1), y = c(0.365,0.25,0.07), hjust = 0.5,size = Fontsize+1)+
  theme(plot.background = element_rect(fill="white", color = NA))+
  draw_plot_label(c("A","B","C","D","E"), c(0.01, 0.45,0.01,0.45,0.8), c(0.995, 0.995, 0.41,0.41,0.41), size = Fontsize) #0.495


ggsave("figure_3.png",f2_cowplot,device="png",width=3200, height=2200,units="px",dpi=300)
ggsave("figure_3.svg",f2_cowplot,device="svg",width=3200, height=2200,units="px")
