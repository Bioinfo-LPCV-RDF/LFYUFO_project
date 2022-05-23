# ## From gene IDs (eg AT1G69120) to gene symbols (eg AP1):
# /home/312.6-Flo_Re/312.6.1-Commun/data/A_thaliana_phytozome_v12/Phytozome/PhytozomeV12/Athaliana/annotation/annotation_AT_TAIR10.csv

library("optparse")
option_list = list(
  make_option(c("-n", "--name"), type="character", default=NULL, 
              help="Set TF and output name", 
              metavar="character"),
  make_option(c("-i", "--in_dir"), type="character",default=NULL, 
              help="Set input directory", 
              metavar="character"),            
  make_option(c("-o", "--out_dir"), type="character",default=NULL, 
              help="Set output directory", 
              metavar="character"), 
  make_option(c("-c", "--control"), type="character",default=NULL, 
              help="Set control set name", 
              metavar="character"), 
  make_option(c("-t", "--treatment"), type="character",default=NULL, 
              help="Set treatment/construct name", 
              metavar="character"), 
  make_option(c("-r", "--reps"), type="character",default=NULL, 
              help="Set number of replicates", 
              metavar="character"),
  make_option(c("-s", "--strategy"), type="numeric",default=NULL, 
              help="Set order of files: treament | control (1) or control | treatment (2)", 
              metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
name=opt$name
dirin=opt$in_dir
dirout=opt$out_dir
control_name=opt$control
treat_name=opt$treatment
nrep=opt$reps
order=as.integer(opt$strategy)


setwd(dirin)

## browseVignettes("package_name") for package information
print("Starting Limma microarray analysis...")
suppressMessages(library(affy))
suppressMessages(library(gcrma))
suppressMessages(library(limma))
suppressMessages(library(ggplot2))
suppressMessages(library(Cairo))
library(calibrate)
library(ggrepel)


## Normalization
## ReadAffy() reads all your cel files
bm=ReadAffy()
bm ## the body map displays the size of the array (feature number), number of samples and number of gen
eset<-gcrma(bm) ## Expression measure, normalised, given in **log2 scale**

setwd(dirout)
write.exprs(eset,file=paste(name,"_expression_estimates.txt",sep="")) ## In this file, you can see the order of your samples, which is useful for your data.cl variable

DATA <- read.table(file=paste(name,"_expression_estimates.txt",sep=""),
                   header=TRUE,sep="\t",row.names = 1)
print(dim(DATA))
colnames(DATA) <- gsub(".CEL","",colnames(DATA))
print(names(DATA))
print("here")
print(order)




# DATA2<-DATA
if (order == 1){
    print("treatment experiment columns followed by control exp columns")
    list<-c(treat_name,control_name)
    Cl<- rep(list, each=nrep)
    print("");print(Cl)
    print(cbind(colnames(DATA),Cl))
    design <- model.matrix(~-1 + factor(Cl))
    colnames(design) <- c(treat_name,control_name)
    print(design)
} else if (order == 2){
    print("control experiment columns followed by treatment exp columns")
    list<-c(control_name,treat_name)
    Cl<- rep(list, each=nrep)
    print("");print(Cl)
    print(cbind(colnames(DATA),Cl))
    design <- model.matrix(~-1 + factor(Cl))
    colnames(design) <- c(control_name,treat_name)
    print(design)
} else if (order == 3){ # 3: control xnrep | treatment xnrep columns
    print("control experiment columns followed by treatment exp columns")
    list<-c(control_name,treat_name)
    Cl<- rep(list, each=nrep)
    print("");print(Cl)
    print(cbind(colnames(DATA),Cl))
    design <- model.matrix(~-1 + factor(Cl))
    colnames(design) <- c(control_name,treat_name)
    print(design)
} else if (order == 4){ # 4: ctrl|treat|ctrl|treat etc, nrep times
    print("control experiment columns followed by treatment exp columns")
    list<-c(treat_name,control_name)
    Cl<- rep(list,nrep)
    print("");print(Cl)
    print(cbind(colnames(DATA),Cl))
    design <- model.matrix(~-1 + factor(Cl))
    colnames(design) <- c(treat_name,control_name)
    print(design)
} else if (order == 5){ # 5: ctrl|treat|ctrl|treat etc, nrep times
    print("control experiment columns followed by treatment exp columns")
#     list<-c(treat_name,control_name)
    Cl<- c("Col","Col","Col","Col","Col","Col","lfy","lfy")
    print("");print(Cl)
    print(cbind(colnames(DATA),Cl))
    design <- model.matrix(~-1 + factor(Cl))
    colnames(design) <- c(treat_name,control_name)
    print(design)
}


fit <- lmFit(DATA, design)

if (order == 1){
    contrast.matrix <- makeContrasts(eval(as.name(control_name)) - eval(as.name(treat_name)), levels = design) # ex for LFYGR: ratio LFYGR/Ler
} else if (order == 2){
    contrast.matrix <- makeContrasts(eval(as.name(treat_name)) - eval(as.name(control_name)), levels = design) # ex for ics1 like in publication
} else if (order == 3){ # 3: treatment/control ratio
    contrast.matrix <- makeContrasts(eval(as.name(control_name)) - eval(as.name(treat_name)), levels = design) # ex for mutants: ratio between wt/lfy
} else if (order == 4){ # 4: LFYGR or 35S/Ler ratio
    contrast.matrix <- makeContrasts(eval(as.name(control_name)) - eval(as.name(treat_name)), levels = design) 
} else if (order == 5){ #5: lfy 7d mutant -> ratio between wt/lfy
    contrast.matrix <- makeContrasts(eval(as.name(treat_name)) - eval(as.name(control_name)), levels = design)
}

fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))

## Volcano plot with unadjusted pvalues
Cairo(width = 700, height = 700, file="ggvolcano_fit.png", type="png", bg = "transparent", canvas = "white", units = "px", dpi = "auto")
volcanoplot(fit2, coef = 1, style = "p-value", xlab = "Log2 Fold Change", ylab = NULL, pch=16, col = rgb(0.6,0.6,0.6,0.6),cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
abline(h=-log10(0.05),col="red")
dev.off()

## MD plot
object<-fit2
Cairo(width = 700, height = 700, file="MDplot.png", type="png", bg = "transparent", canvas = "white", units = "px", dpi = "auto")
plotMD(object, column = ncol(object), coef = NULL, xlab = "Average log-expression",
       ylab = "log-fold-change", main = "MD plot of fit",
       status=object$genes$Status, zero.weights = FALSE,pch=16,col =rgb(0,0,0,0.6),cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
abline(h=0,col="red",lty=3)
dev.off()

## MDS plot
dim.plot = c(1,2)
Cairo(width = 700, height = 700, file="MDSplot.png", type="png", bg = "transparent", canvas = "white", units = "px", dpi = "auto")
plotMDS(DATA, top = 500, labels = NULL, pch = NULL, cex = 1,
        dim.plot = c(1,2), ndim = max(dim.plot), gene.selection = "pairwise",
        xlab = NULL, ylab = NULL, plot = TRUE,cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)
dev.off()



fit2_EZH2 <- fit2 # pour l'utiliser plus tard 
res1 <- topTable(fit2, coef = 1, adjust = "BH", number = nrow(DATA))  # extraire resultats de la 1ere comparaison -- already adjusted!
print(head(res1))
ath1_probe<-rownames(res1)
rownames(res1)<-NULL
res1<-cbind(ath1_probe,res1)
head(res1)
write.table(res1,file="results1.csv",sep=';',row.names=F)


### add ATXG codes ###
gene_descr_file<-"/home/312.6-Flo_Re/312.6.1-Commun/LFY/LFY_targets/data/short_description.csv"
gene_descr <- read.csv(gene_descr_file, header =TRUE, sep =";",quote ='"')
print("gene description file successfully uploaded")

names(gene_descr)[1]<-"ath1_probe"
gene_descr<-gene_descr[,c(1,2)]
print(str(gene_descr))
res_limma<-read.csv("results1.csv", header =TRUE, sep =";",quote ='"')


df<-merge(res_limma,gene_descr,by=("ath1_probe"),all=T)
df$Locus.Identifier<-as.character(df$Locus.Identifier)
df$ath1_probe<-as.character(df$ath1_probe)


df<-na.omit(df)
print("");print(dim(df))


#### Volcano plots ####
df$diffexpressed <- "NO"
df$diffexpressed[df$Locus.Identifier == "AT5G03840"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT5G03790"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT3G61250"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT5G23000"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT1G69120"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT3G54340"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT5G24910"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT4G18960"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT5G15800"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT3G02310"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT1G24260"]<- "poscon"
df$diffexpressed[df$Locus.Identifier == "AT5G20240"]<- "poscon"

df$diffexpressed[df$Locus.Identifier == "AT5G61850"] <- "LFY"

print("1")
df$delabel <- NA
df$delabel[df$diffexpressed == "LFY"] <- "LFY"

df$delabel[df$Locus.Identifier == "AT5G03840"]<- "TFL1"
df$delabel[df$Locus.Identifier == "AT5G03790"]<- "LMI1"
df$delabel[df$Locus.Identifier == "AT3G61250"]<- "LMI2"
df$delabel[df$Locus.Identifier == "AT5G23000"]<- "RAX1"
df$delabel[df$Locus.Identifier == "AT1G69120"]<- "AP1"
df$delabel[df$Locus.Identifier == "AT3G54340"]<- "AP3"
df$delabel[df$Locus.Identifier == "AT5G24910"]<- "ELA1"
df$delabel[df$Locus.Identifier == "AT4G18960"]<- "AG"
df$delabel[df$Locus.Identifier == "AT5G15800"]<- "SEP1"
df$delabel[df$Locus.Identifier == "AT3G02310"]<- "SEP2"
df$delabel[df$Locus.Identifier == "AT1G24260"]<- "SEP3"
df$delabel[df$Locus.Identifier == "AT5G20240"]<- "PI"


print("1")
print("and now create volcano plot!")
volcanodata<-as.data.frame(cbind(df$ath1_probe,df$Locus.Identifier,df$logFC,df$adj.P.Val,df$B))

names(volcanodata)<-c("ath1_probe","Locus.Identifier","logFC","adj_pvalue","B")

volcanodata$logFC<-as.numeric(as.character(volcanodata$logFC))
volcanodata$adj_pvalue<-as.numeric(as.character(volcanodata$adj_pvalue))
volcanodata$B<-as.numeric(as.character(volcanodata$B))

write.table(volcanodata,file=paste(name,"_volcanodata.csv",sep=""),sep=';',row.names=F)


DEG<-subset(volcanodata,(logFC>0.5 & adj_pvalue<0.05) | (logFC<(-0.5) & adj_pvalue<0.05))
write.table(DEG,file=paste(name,"_DEG.csv",sep=""),sep=';',row.names=F)

DEG_FC1<-subset(volcanodata,(logFC>1 & adj_pvalue<0.05) | (logFC<(-1) & adj_pvalue<0.05))
write.table(DEG_FC1,file=paste(name,"_DEG_FC1-1.csv",sep=""),sep=';',row.names=F)


#### Plot volcanos based on df data (not volcano plot data) ####
mycolors <- c("yellow", "red", "darkgrey")
names(mycolors) <- c("LFY", "poscon", "NO")

Cairo(width = 700, height = 700, file="ggvolcano.png", type="png", bg = "transparent", canvas = "white", units = "px", dpi = "auto")

ggplot(data=df, aes(x=logFC, y=-log10(adj.P.Val), 
           col=diffexpressed,label=delabel))+
    geom_point(size=ifelse(df$diffexpressed == 'NO',1,2), 
               alpha=ifelse(df$diffexpressed == 'NO',0.3,1))+
    ggtitle(name)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 20,face = "bold"),
          legend.position = "none",
          axis.title=element_text(size=22),
          axis.text=element_text(size=20))+
    scale_colour_manual(values = mycolors) +
    geom_hline(yintercept=-log10(0.05), col="red",
               size=0.1)+
    geom_vline(xintercept =c(-1,-0.5,0.5,1), col = c("green","blue","blue","green"), size=0.2)+
    geom_text_repel(colour = "black",size=8)+
    xlim(-6,6)+
    ylim(0,6)

dev.off()










