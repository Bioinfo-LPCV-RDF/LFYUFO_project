## script to plot PI promoter sequence sites of m and dLUBS
library(ggplot2)
library(svglite)
args <- commandArgs(TRUE)
wd=args[1]
name=args[2]
print(paste0("directory where data is stored: ",wd))
setwd(wd)

# setwd("//oslo/share_nfs/312-Projets_PCV/312.6-Flo_Re/312.6.1-Commun/LFY/LFY-UFO/results/comparisons/LFYamp_LFYUFO_cf/RBE")
# setwd(wd)

name="RBE"

### dLUBS
dLUBS<-read.table(paste0(name,"_dLUBS_light.scores"),sep="\t",header=F,dec=".")
# head(dLUBS)
names(dLUBS)<-c("pos","strand","score")


dLUBS_subset<-subset(dLUBS,score>=-13)## 2 | V2 < 4)
print(dLUBS_subset)

vlines <- function(x, y, y.orig,y.max,col){
  plot(range(x), range(c(y, y.orig)), type='n',
       ylim = c(y.orig,y.max),xlim = c(0,1562),#c(-975,0),
       yaxt="n",xaxt="n",
       axes = FALSE,
       ylab="score",xlab="",cex.lab=1.2)
  for (i in seq(along=x))lines(x[c(i,i)], c(y.orig, y[i]),col = col,lwd=6,lend=1)
  axis(side=2, at=c(y.max,y.orig), labels = FALSE)
  text(par("usr")[1], c(y.max,y.orig),  
       labels = c(y.max,y.orig), pos = 2, xpd = T)
  # custom x axis
  xtick<-c(0,1562) #seq(0,-975, by=-500)
  axis(side=1, at=xtick, labels = FALSE)
  text(x=xtick,  par("usr")[3], 
       labels = c("ATG","1562"),
       pos = 1, xpd = TRUE)
  
}


png(paste0(name,"_prom_scores_dLUBS.png"),width = 2000,height = 800,units = "px",res = 300)
vlines(x=dLUBS_subset$pos, y=dLUBS_subset$score, y.orig=-13,y.max=-7,"red")
dev.off()

svg(paste0(name,"_prom_scores_dLUBS.svg"),width = 9,height = 2.5)#,units = "px")#,res = 300)
vlines(x=dLUBS_subset$pos, y=dLUBS_subset$score, y.orig=-13,y.max=-7,"red")
dev.off()






### mLUBS
mLUBS<-read.table(paste0(name,"_mLUBS_light.scores"),sep="\t",header=F,dec=".")
# head(mLUBS)
names(mLUBS)<-c("pos","strand","score")


mLUBS_subset<-subset(mLUBS,score>=-9)## 2 | V2 < 4)
print(mLUBS_subset)

png(paste0(name,"_prom_scores_mLUBS.png"),width = 2000,height = 800,units = "px",res = 300)
vlines(x=mLUBS_subset$pos, y=mLUBS_subset$score, y.orig=-9,y.max=-4,"blue")
dev.off()

svg(paste0(name,"_prom_scores_mLUBS.svg"),width = 9,height = 2.5)#,units = "px")#,res = 300)
vlines(x=mLUBS_subset$pos, y=mLUBS_subset$score, y.orig=-9,y.max=-4,"blue")
dev.off()








## LFY
LFY<-read.table(paste0(name,"_LFY_light.scores"),sep="\t",header=F,dec=".")
# head(LFY)
names(LFY)<-c("pos","strand","score")


LFY_subset<-subset(LFY,score>=-21)

png(paste0(name,"_prom_scores_LFY.png"),width = 2000,height = 800,units = "px",res = 300)
vlines(x=LFY_subset$pos, y=LFY_subset$score, y.orig=-21,y.max=-12,"forestgreen")
dev.off()

svg(paste0(name,"_prom_scores_LFY.svg"),width = 9,height = 2.5)#,units = "px")#,res = 300)
vlines(x=LFY_subset$pos, y=LFY_subset$score, y.orig=-21,y.max=-12,"forestgreen")
dev.off()




quit()










quit()

plot(mLUBS_subset$pos,mLUBS_subset$score, bty="n",type = "h",
     col="blue",lwd=2,lend=1,
     ylim=c(-4,-9),
     xlim=c(-1500,0),
     yaxt="n",xaxt="n",
     ylab="score",xlab="",cex.lab=1.2)
# custom y axis
axis(side=2, at=c(-4,-9), labels = FALSE)
text(par("usr")[1], c(-4,-9),  
     labels = c(-4,-9), pos = 2, xpd = T)
# custom x axis
xtick<-seq(-1500, 0, by=500)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3], 
     labels = c(-1500,-1000,-500,"ATG"), pos = 1, xpd = TRUE)






