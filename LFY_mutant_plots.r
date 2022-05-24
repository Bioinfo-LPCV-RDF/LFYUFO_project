## script for all plots of LFY, LFYUFO, LFYm and LFYmUFO
library(ggplot2)
library(svglite)
library(reshape2)



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
print(wilcox.test(LFYUFO_spe_peaks$CFC,LFYUFO_spe_peaks$LFYmLFY))


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


ggsave(file="violinplot_top20perc_CFC_pub.png",device="png",plot=violin_ratios_pub,units="px", width=1800, height=1500)
ggsave(file="violinplot_top20perc_CFC_pub.svg",device="svg",plot=violin_ratios_pub,units="px", width=1800, height=1500)




quit()


