## script for all plots of LFY, LFYUFO, LFYm and LFYmUFO
library(ggplot2)
library(svglite)
library(reshape2)


just a change 

args <- commandArgs(TRUE)
wd_comp1=args[1]
wd_comp2=args[2]
wd_comp3=args[3]

print(paste0("setting wd to: ",wd_comp1))


## Comparison 1: LFYm vs LFYmUFO
setwd(wd_comp1)
table<-read.table("peaks_perSample_rpkminPeaks.txt",sep="\t",header=T,dec=".")
head(table)

table$LFYm<-rowMeans(table[,c('LFYm_a', 'LFYm_b')], na.rm=TRUE)
table$LFYmUFO<-rowMeans(table[,c('LFYmUFO_cfa', 'LFYmUFO_cfb')], na.rm=TRUE)

table$CFCm<-(table$LFYmUFO/table$LFYm)




LFYmvsLFYmUFO_gradient=c("#4cccff", "#a5a4a5", "#c39788", "#d29079", "#f0835b", "#ff7c4c") # simplified





Fontsize=10
limitsplot=c(min(c(table$LFYmUFO,table$LFYm)),max(c(table$LFYmUFO,table$LFYm)))

## graph 1: LFYm vs LFYmUFO, colored by LFYmUFO/LFYm ratio
p_LFYm_LFYmUFO<-ggplot(data=table,aes(x=LFYmUFO, y=LFYm))+
  scale_y_log10(limits=limitsplot)+
  scale_x_log10(limits=limitsplot)+
  geom_point(alpha=0.3,aes(color=CFCm))+
  geom_abline(intercept = 0,linetype="dashed",size=0.1) +
  scale_colour_gradientn(colours=LFYmvsLFYmUFO_gradient,name=expression(paste("CFC = ",frac(paste(LFY[K249R],"+UFO"),LFY[K249R]))))+
  theme(text = element_text(size = 7))+
  xlab(expression(paste(LFY[K249R],"+UFO (",Log[10]," RPKM)")))+
  ylab(expression(paste(LFY[K249R]," (",Log[10]," RPKM)")))+
  theme_classic()+coord_fixed()+
  theme(axis.title=element_text(size=Fontsize+4),
        axis.text.x=element_text(size=Fontsize+2),
        axis.text.y=element_text(size=Fontsize+2),
        legend.text=element_text(size=Fontsize+3),
        legend.title=element_text(size=Fontsize+3),
        legend.box = "horizontal")


ggsave(file="CFC_LFYm_LFYmUFO.svg",device="svg",plot=p_LFYm_LFYmUFO,units="px", width=2000, height=1500)
ggsave(file="CFC_LFYm_LFYmUFO.png",device="png",plot=p_LFYm_LFYmUFO,units="px", width=2000, height=1500)





## Comparison 2: LFYamp vs LFYm
print(paste0("setting wd to: ",wd_comp2))
setwd(wd_comp2)
table<-read.table("peaks_perSample_rpkminPeaks.txt",sep="\t",header=T,dec=".")
head(table)

table$LFYm<-rowMeans(table[,c('LFYm_a', 'LFYm_b')], na.rm=TRUE)
table$LFYamp<-rowMeans(table[,c('LFYamp1', 'LFYamp2', 'LFYamp3')], na.rm=TRUE)

table$LFYmLFY<-(table$LFYm/table$LFYamp)

Fontsize<-12
limitsplot<-c(min(c(table$LFYm,table$LFYamp)),max(c(table$LFYm,table$LFYamp)))


## graph 2: covplot LFY vs LFYm
p_LFYm_LFY<-ggplot(data=table,aes(x=LFYamp, y=LFYm))+
  scale_y_log10(limits=limitsplot)+
  scale_x_log10(limits=limitsplot)+
  geom_point(alpha=0.3,color="#931F1D")+
  geom_abline(intercept = 0,linetype="longdash",size=0.1) +
  theme(text = element_text(size = 7))+
  xlab(expression(paste("LFY (",Log[10]," RPKM)")))+
  ylab(expression(paste(LFY[K249R]," (",Log[10]," RPKM)")))+
  theme_classic()+coord_fixed()+
  theme(axis.title=element_text(size=Fontsize+4),
        axis.text.x=element_text(size=Fontsize+2),
        axis.text.y=element_text(size=Fontsize+2),
        legend.text=element_text(size=Fontsize+3),
        legend.title=element_text(size=Fontsize+3),
        legend.box = "horizontal")


ggsave(file="CFC_LFY_LFYm.svg",device="svg",plot=p_LFYm_LFY,units="px", width=1500, height=1500)
ggsave(file="CFC_LFY_LFYm.png",device="png",plot=p_LFYm_LFY,units="px", width=1500, height=1500)




### look at effect of mutation on top 20% CFC peaks
print(paste0("setting wd to: ",wd_comp3))
setwd(wd_comp3)
table<-read.table("peaks_perSample_rpkminPeaks.txt",sep="\t",header=T,dec=".")
print(dim(table))

## define columns
table$LFYamp<-rowMeans(table[,c('LFYamp1', 'LFYamp2', 'LFYamp3')], na.rm=TRUE)
table$LFYUFO<-rowMeans(table[,c('LFYUFO_cfa', 'LFYUFO_cfb', 'LFYUFO_cfc')], na.rm=TRUE)
table$LFYm<-rowMeans(table[,c('LFYm_a', 'LFYm_b')], na.rm=TRUE)
table$LFYmUFO<-rowMeans(table[,c('LFYmUFO_cfa', 'LFYmUFO_cfb')], na.rm=TRUE)

## define ratios
table$CFC<-(table$LFYUFO/table$LFYamp)
table$CFCm<-(table$LFYmUFO/table$LFYm)
table$LFYmLFY<-(table$LFYm/table$LFYamp)
print(dim(table))

CFC_sorted_table<-table[order(-table$CFC),]
twentyperc<-round(dim(CFC_sorted_table)[1]*0.2,digits=0)
LFYUFO_spe_peaks<-head(CFC_sorted_table, n=twentyperc)


LFYUFO_spe_peaks<-LFYUFO_spe_peaks[,c(14:20)] # only keep last six columns
print(head(LFYUFO_spe_peaks))



## test significance of difference in median
print(wilcox.test(LFYUFO_spe_peaks$CFC,LFYUFO_spe_peaks$CFCm))
print(wilcox.test(LFYUFO_spe_peaks$CFCm,LFYUFO_spe_peaks$LFYmLFY))
wilcox.test(LFYUFO_spe_peaks$CFC,LFYUFO_spe_peaks$LFY_LFYm)


### figure for publication 
melted_top20_CFC = melt(LFYUFO_spe_peaks, id.vars = c("CFCm","LFYmLFY"), 
                        measure.vars = c("CFC", "CFCm","LFYmLFY"))



violin_ratios_pub<-ggplot(data=melted_top20_CFC,aes(x=as.factor(variable), 
                                                    y=value,
                                                    fill=as.factor(variable)))+
  scale_y_log10()+
  geom_violin(trim=FALSE)+
  geom_boxplot(width=.1, position = position_dodge(0.9),col = "grey20",fill="white")+
  ylab(expression(paste("ratio value (",Log[10],")")))+
  theme_classic()+
  theme(axis.title.y=element_text(size=Fontsize+4),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=Fontsize+2),
        axis.text.y=element_text(size=Fontsize+2),
        legend.position="none")+
  scale_fill_manual(values=c("#fcf6f6", "#ebbab9", "#da7e7c"),
                    name = "ratio", labels = c(expression(paste("CFC = ",frac("LFY+UFO","LFY"))),
                                               expression(paste("CFC = ",frac(paste(LFY[K249R],"+UFO"),LFY[K249R]))),
                                               expression(paste(LFY[K249R],"/LFY"))))+
  scale_x_discrete(labels=c(expression(paste(frac("LFY+UFO","LFY"))),
                            expression(paste(frac(paste(LFY[K249R],"+UFO"),LFY[K249R]))),
                            expression(paste(frac(LFY[K249R],"LFY")))))

# violin_ratios_pub

# "#e9f5f8", "#add8e6", "#75a9d7"
# "#F4EEFF", "#DCD6F7", "#A6B1E1"
# "#BBDEF0", "#6B2737", "#E08E45"
# "#FFBF46","#4B3B47","#6A6262"
# "#59C3C3", "#52489C", "#FFBF46"
# "#F9A12EFF","#FC766AFF","#9B4A97FF"

ggsave(file="violinplot_top20perc_CFC_pub.png",device="png",plot=violin_ratios_pub,units="px", width=1800, height=1500)
ggsave(file="violinplot_top20perc_CFC_pub.svg",device="svg",plot=violin_ratios_pub,units="px", width=1800, height=1500)




quit()


















## same on LFY-spe sites, so bottom 20% CFC
twentyperc<-round(dim(CFC_sorted_table)[1]*0.20,digits=0)
LFY_spe_20<-tail(CFC_sorted_table, n=twentyperc)


LFY_spe_20<-LFY_spe_20[,c(14:22)] ## only keep last six columns


## melted top 20%
melted_bot20_cov = melt(LFY_spe_20, id.vars = c("CFC","CFC_dec","CFCm","LFY_LFYm","LFYmLFY"),
                        measure.vars = c("LFYUFO", "LFYmUFO","LFYamp","LFYm"))

head(melted_bot20_cov)

violin_cov_bot20<-ggplot(data=melted_bot20_cov,aes(x=as.factor(variable), y=value))+#,group=interaction(CFC_dec,variable)))+
  scale_y_log10()+#limits=limitsplot)+
  geom_violin(aes(fill=variable),trim=FALSE)+#,position = "dodge")+
  geom_boxplot(width=.1, position = position_dodge(0.9),col = "grey20")+ ## outlier.colour=NA,
  ylab("log10(RPKM) on bottom 20% peaks")+xlab("dataset")+
  theme_classic()


melted_bot20_CFC = melt(LFY_spe_20, id.vars = c("CFC_dec","CFCm","LFYmLFY"),
                        measure.vars = c("CFC", "CFCm","LFYmLFY","LFY_LFYm"))#,"LFYamp","LFYm"))
head(melted_bot20_CFC)


violin_ratios_bot20<-ggplot(data=melted_bot20_CFC,aes(x=as.factor(variable), y=value))+#,group=interaction(CFC_dec,variable)))+
  scale_y_log10()+#limits=limitsplot)+
  geom_violin(aes(fill=variable),trim=FALSE)+
  geom_boxplot(width=.1, position = position_dodge(0.9),col = "grey20")+ ## outlier.colour=NA,
  xlab("dataset")+ylab(expression(paste("ratio value (",Log[10],")")))+
  # scale_fill_manual(values=c("blue","orange","pink"))+
  theme_classic()+
  theme(axis.title=element_text(size=Fontsize+4),
        axis.text.x=element_text(size=Fontsize+2),
        axis.text.y=element_text(size=Fontsize+2),
        legend.text=element_text(size=Fontsize+3),
        legend.title=element_text(size=Fontsize+3),
        legend.box = "horizontal")





















## same with 15% top CFC
fifteenperc<-round(dim(CFC_sorted_table)[1]*0.15,digits=0)
LFYUFO_spe_15<-head(CFC_sorted_table, n=fifteenperc)


LFYUFO_spe_15<-LFYUFO_spe_15[,c(14:22)] ## only keep last six columns


## melted top 20%
melted_top15_cov = melt(LFYUFO_spe_15, id.vars = c("CFC","CFC_dec","CFCm","LFY_LFYm","LFYmLFY"),
                        measure.vars = c("LFYUFO", "LFYmUFO","LFYamp","LFYm"))

head(melted_top15_cov)

violin_cov_15<-ggplot(data=melted_top15_cov,aes(x=as.factor(variable), y=value))+#,group=interaction(CFC_dec,variable)))+
  scale_y_log10()+#limits=limitsplot)+
  geom_violin(aes(fill=variable),trim=FALSE)+#,position = "dodge")+
  geom_boxplot(width=.1, position = position_dodge(0.9),col = "grey20")+ ## outlier.colour=NA,
  ylab("log10(RPKM) on top 20% peaks")+xlab("dataset")+
  theme_classic()


melted_top15_CFC = melt(LFYUFO_spe_15, id.vars = c("CFC_dec","CFCm","LFYmLFY"),
                        measure.vars = c("CFC", "CFCm","LFYmLFY","LFY_LFYm"))#,"LFYamp","LFYm"))
head(melted_top15_CFC)


violin_ratios_15<-ggplot(data=melted_top15_CFC,aes(x=as.factor(variable), y=value))+#,group=interaction(CFC_dec,variable)))+
  scale_y_log10()+#limits=limitsplot)+
  geom_violin(aes(fill=variable),trim=FALSE)+
  geom_boxplot(width=.1, position = position_dodge(0.9),col = "grey20")+ ## outlier.colour=NA,
  xlab("dataset")+ylab(expression(paste("ratio value (",Log[10],")")))+
  # scale_fill_manual(values=c("blue","orange","pink"))+
  theme_classic()+
  theme(axis.title=element_text(size=Fontsize+4),
        axis.text.x=element_text(size=Fontsize+2),
        axis.text.y=element_text(size=Fontsize+2),
        legend.text=element_text(size=Fontsize+3),
        legend.title=element_text(size=Fontsize+3),
        legend.box = "horizontal")






##### rep/rep plots LFYm and LFYmUFO #####
limitsplot_LFYm=c(min(table$LFYm_a,table$LFYm_b),max(table$LFYm_a,table$LFYm_b))
LFYm<-ggplot(table,aes(x=LFYm_a,y=LFYm_b))+
  scale_y_log10(limits=limitsplot_LFYm)+
  scale_x_log10(limits=limitsplot_LFYm)+
  geom_point(alpha=0.3,col="grey45")+
  geom_abline(intercept = 0,linetype="longdash",size=0.1)+
  xlab(expression(paste(LFY[K249R]," rep1 (",Log[10]," RPKM)")))+
  ylab(expression(paste(LFY[K249R]," rep2 (",Log[10]," RPKM)")))+
  theme_classic()+coord_fixed()+
  theme(axis.title=element_text(size=11),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        plot.margin = unit(c(1,15,1,1),"pt"))
ggsave(filename ="LFYm_RIP_rep_plot.svg",device="svg",plot=LFYm,width=1000,height=1000,unit="px")
ggsave(filename ="LFYm_RIP_rep_plot.png",device="png",plot=LFYm,width=1000,height=1000,unit="px")


limitsplot_LFYmUFO=c(0.1,max(table$LFYmUFO_cfa,table$LFYmUFO_cfb))
LFYmUFO<-ggplot(table,aes(x=LFYmUFO_cfa,y=LFYmUFO_cfb))+
  scale_y_log10(limits=limitsplot_LFYmUFO)+
  scale_x_log10(limits=limitsplot_LFYmUFO)+
  geom_point(alpha=0.3,col="grey45")+
  geom_abline(intercept = 0,linetype="longdash",size=0.1)+
  xlab(expression(paste(LFY[K249R],"+UFO rep1 (",Log[10]," RPKM)")))+
  ylab(expression(paste(LFY[K249R],"+UFO rep2 (",Log[10]," RPKM)")))+
  theme_classic()+coord_fixed()+
  theme(axis.title=element_text(size=11),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        plot.margin = unit(c(1,15,1,1),"pt"))
ggsave(filename ="LFYmUFO_RIP_rep_plot.svg",device="svg",plot=LFYmUFO,width=1000,height=1000,unit="px")
ggsave(filename ="LFYmUFO_RIP_rep_plot.png",device="png",plot=LFYmUFO,width=1000,height=1000,unit="px")






