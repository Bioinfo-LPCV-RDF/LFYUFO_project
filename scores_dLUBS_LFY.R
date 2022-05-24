## script to plot figure S##: comparison of LFYUFO vs LFY scores on LFYUFO sites
library(ggplot2)#,lib.loc="/home/312.3-StrucDev/312.3.1-Commun/R-3.5.0/")
options(bitmapType='cairo') #keep this line
library(reshape2)
library(cowplot)


args <- commandArgs(TRUE)
print(args)
wd=args[1]

setwd(wd)
print("Creating boxplot for LFY quality assessment")

tableName<-"table_all_scores.tsv"

table <- read.table(file=tableName, header=TRUE, sep = "", dec=".", stringsAsFactors=TRUE)

table$Decile1 <- with(table, factor(
  findInterval( CFC, c(-Inf,
                       quantile(CFC, probs=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)), Inf) ), 
  labels=c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")))

table$decile <- with(table, cut(CFC, 
                                breaks=quantile(CFC, probs=seq(0,1, by=0.1)), 
                                include.lowest=TRUE))


## same thing but only top 20% CFC, boxplot (or violin?)
twentyperc<-round(((dim(table)[1])*0.2),digits=0)

## sort table by decreasing CFC value
table_CFC_ordered<-table[order(table$CFC,decreasing=TRUE),]

## select top 20% CFC peaks
top20CFC<-head(table_CFC_ordered,n=twentyperc)
top20CFC<-top20CFC[,c(4,7:10)]

melted<-melt(top20CFC, id=c("file","CFC"))


## select lowest 20% of CFC
low20CFC<-tail(table_CFC_ordered,n=twentyperc)
low20CFC<-low20CFC[,c(4,7:10)]

melted_lowCFC<-melt(low20CFC, id=c("file","CFC"))


## merge top and bottom CFC to have: 
## LFY_on_dLUBS values for top CFC and
## LFY values for bottom CFC
df1<-subset(melted,variable=="LFY_on_dLUBS") # df with only LFY_on_dLUBS values
df2<-subset(melted_lowCFC,variable=="LFY") # df with only LFY values

merged<-rbind(df1,df2)
head(merged)

boxplot_mixed<-ggplot(merged,aes(x=variable,y=value,fill=variable))+
  geom_boxplot()+
  ylim(-52,-4)+
  ylab("score")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        axis.title.x = element_blank())

ggsave("boxplots_dLUBS_vs_LFY.png",boxplot_mixed,width=1500,height = 1200,unit="px")
ggsave("boxplots_dLUBS_vs_LFY.svg",boxplot_mixed,width=1500,height = 1200,unit="px")

quit()

wilcox.test(df1$value,df2$value)


