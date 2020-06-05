library(data.table)
library(metagenomeSeq)
library(phyloseq)

# Run with --help flag for help.
# Written 02/05/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="basepath/tables/Taxa_merged_species.txt",
              help="Input file", metavar="character"),
  make_option(c("-A", "--analysis"), type="character", default="Peressotti",
              help="Output file", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="basepath/plots/",
              help="Output file", metavar="character"),
  make_option(c("-R", "--rank"), type="character", default="Family",
              help="Rank to be used for analysisi", metavar="character")
 # make_option(c("-T", "--treated_prefix"), type="character", default="experimental_",
              # help="Prefix identifying treated samples", metavar="character"),
  # make_option(c("-C", "--control_prefix"), type="character", default="control_",
              # help="Prefix identifying control samples", metavar="character"),
  # make_option(c("-R", "--raw_counts"), type="character", default=NULL,
              # help="raw (total) read counts for this starting file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No infile file specified with '-I' flag.")
} else {  cat ("infile is ", opt$infile, "\n")
  infile <- opt$infile  
  }

  if (is.null(opt$outdir)) {
  stop("WARNING: No output directory specified with '-O' flag.")
} else {  cat ("Output directory is ", opt$outdir, "\n")
  outdir <- opt$outdir  
  }

  if (is.null(opt$rank)) {
  stop("WARNING: No rank specified with '-R' flag.")
} else {  cat ("Rank is ", opt$rank, "\n")
  rank <- opt$rank  
  }

tri.to.squ<-function(x)
{
rn<-row.names(x)
cn<-colnames(x)
an<-unique(c(cn,rn))
myval<-x[!is.na(x)]
mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
for(ext in 1:length(cn))
{
 for(int in 1:length(rn))
 {
 if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
 mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 }
  
}
return(mymat)
}




plot.diversity<-function(infile,outdir,rank)
{
library("data.table")
library("vegan")
library(multcompView)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
mydata<-fread(paste("basepath/tables/16s_merged_count_by_",rank,".txt",sep=""),data.table=F)
browser()
#Set up a metadata on the fly
if(analysis=="Basic")
{
Group<-c(35,35,35,18,18,18,9,9,9,4,4,4,0,0,0,"FTQ","FD")
sample<-seq(1,17)
metadata<-data.frame(sample,Group,check.names=F)
#Set up the comparison (here 35psu vs 0)
tocompare<-c(35,0)
samplecomp<-metadata$sample[metadata$Group%in%tocompare]
metadata<-metadata[metadata$Group%in%tocompare,]
}


if(analysis=="Peressotti")
{
Group<-PerCla<-c("A","A","A","A","A","A",rep("B",9),"FTQ","FD")
sample<-seq(1,17)
metadata<-data.frame(sample,Group,check.names=F)
#Set up the comparison (here 35psu vs 0)
tocompare<-c("A","B")
samplecomp<-metadata$sample[metadata$Group%in%tocompare]
metadata<-metadata[metadata$Group%in%tocompare,]
}

row.names(mydata)<-mydata$V1
mydata$Tot<-mydata$V1<-NULL
mydata<-mydata[,seq(1:15)]
mydata<-t(mydata)

sim.div<-diversity(mydata,"simpson")
shan.div<-diversity(mydata,"shannon")
species.richness<-estimateR(mydata)
pino<-data.frame(Names=names(sim.div),Simpson=round(sim.div,4),Shannon=round(shan.div,4),Species_Richness=species.richness[1,])
pino$salinity<-c(35,35,35,18,18,18,9,9,9,4,4,4,0,0,0)
write.table(pino,paste("basepath/tables/alpha_diversity_",rank,".txt",sep=""),sep="\t",quote=F,row.names=F)
#Perform and write spearman correlation test: we perfomr on all the values, and collapsing (averaging) biological replicates, to see if results change (they don't!)
fullres<-data.frame(Test=c("Full_simpson","Full_shannon","Full_chao","Mean_simpson","Mean_shannon","Mean_chao"),rho=0,pvalue=1)

si<-cor.test(pino$Simpson,pino$salinity,method="spearman")
sh<-cor.test(pino$Shannon,pino$salinity,method="spearman")
sr<-cor.test(pino$Species_Richness,pino$salinity,method="spearman")
cc<-aggregate(pino[,2:ncol(pino)],by=list(pino$salinity),FUN="mean")
msi<-cor.test(cc$Simpson,cc$Group.1,method="spearman")
msh<-cor.test(cc$Shannon,cc$Group.1,method="spearman")
msr<-cor.test(cc$Species_Richness,cc$Group.1,method="spearman")

fullres[1,2]<-si$estimate
fullres[1,3]<-si$p.value
fullres[2,2]<-sh$estimate
fullres[2,3]<-sh$p.value
fullres[3,2]<-sr$estimate
fullres[3,3]<-sr$p.value
fullres[4,2]<-msi$estimate
fullres[4,3]<-msi$p.value
fullres[5,2]<-msh$estimate
fullres[5,3]<-msh$p.value
fullres[6,2]<-msr$estimate
fullres[6,3]<-msr$p.value

write.table(fullres,paste("basepath/tables/corr_alpha_diversity_salinity_",rank,".txt",sep=""),sep="\t",quote=F,row.names=F)


browser()
pdf(paste(outdir,"alpha_diversity.pdf",sep=""),width=14,height=6)
par(mfrow=c(1,3))
barplot(sim.div,main="Simpson diversity")
barplot(shan.div,main="Shannon's Index")
barplot(species.richness[1,],main="Species richness")
dev.off()
#tiff(paste(outdir,"alpha_diversity.tiff",sep=""),width=14,height=6)
#par(mfrow=c(1,3))
#barplot(sim.div,main="Simpson diversity")
#barplot(shan.div,main="Shannon's Index")
#barplot(species.richness[1,],main="Species richness")
#dev.off()
pdf(paste(outdir,"simpson_div.pdf",sep=""))
barplot(sim.div)
dev.off()
pdf(paste(outdir,"shannon_div.pdf",sep=""))
barplot(shan.div)
dev.off()
pdf(paste(outdir,"species_richness.pdf",sep=""))
barplot(species.richness[1,])
dev.off()
}




fake<-function()
{
fungi<-aggregate(fungi[,!names(fungi)%in%"binomia"],by=list(fungi$binomia),FUN="sum")
oomyc<-aggregate(oomyc[,!names(oomyc)%in%"binomia"],by=list(oomyc$binomia),FUN="sum")
row.names(fungi)<-fungi$"Group.1"
row.names(oomyc)<-oomyc$"Group.1"
fungi$"Group.1"<-NULL
oomyc$"Group.1"<-NULL
fungi<-t(fungi)
oomyc<-t(oomyc)
fungiplants<-fungi[1:24,]
fungisoil<-fungi[25:48,]
oomycplants<-oomyc[1:24,]
oomycsoil<-oomyc[25:48,]

#Compute simpson diversity
f_p_simpson <- diversity(fungiplants, "simpson")
o_p_simpson <- diversity(oomycplants, "simpson")
f_s_simpson <- diversity(fungisoil, "simpson")
o_s_simpson <- diversity(oomycsoil, "simpson")
#Pairwise wilcoxon test
simdata<-data.frame(Substrate=c(rep("Plants\nFungi",24),rep("Plants\nOomycetes",24),rep("Soil\nFungi",24),rep("Soil\nOomycetes",24)),
Diversity=c(f_p_simpson,o_p_simpson,f_s_simpson,o_s_simpson))
psimps<-pairwise.wilcox.test(simdata$Diversity, simdata$Substrate, p.adjust.method = "none", paired = FALSE)
simpsmat<-tri.to.squ(psimps$p.value)
simletters<-multcompLetters(simpsmat,compare="<=",threshold=0.05,Letters=letters)$Letters

#Compute shannon diversity
f_p_shannon <- diversity(fungiplants, "shannon")
o_p_shannon <- diversity(oomycplants, "shannon")
f_s_shannon <- diversity(fungisoil, "shannon")
o_s_shannon <- diversity(oomycsoil, "shannon")
#Pairwise wilcoxon test
shadata<-data.frame(Substrate=c(rep("Plants\nFungi",24),rep("Plants\nOomycetes",24),rep("Soil\nFungi",24),rep("Soil\nOomycetes",24)),
Diversity=c(f_p_shannon,o_p_shannon,f_s_shannon,o_s_shannon))
pshan<-pairwise.wilcox.test(shadata$Diversity, shadata$Substrate, p.adjust.method = "none", paired = FALSE)
shamat<-tri.to.squ(pshan$p.value)
shaletters<-multcompLetters(shamat,compare="<=",threshold=0.05,Letters=letters)$Letters



mycol<-brewer.pal(n=4,"Set1")



#Plot alpha diversity stratifying by substrate (plants/soil)
#png(paste(outdir,"diversity_by_substrate.pdf",sep=""),width=6,height=6,units="cm",res=600,type="cairo")
pdf(paste(outdir,"diversity_by_substrate.pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2))
  pino<-ggplot(simdata, aes(x=Substrate, y=Diversity)) +
  xlab("Sample")+ylab("Simpson Diversity")+ 
  geom_boxplot(outlier.shape=NA,fill=mycol) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0)) +
  geom_text(aes(x=1, y=max(f_p_simpson)+0.05*max(f_p_simpson),label=simletters[1])) +
  geom_text(aes(x=2, y=max(o_p_simpson)+0.05*max(o_p_simpson),label=simletters[2])) +
  geom_text(aes(x=3, y=max(f_s_simpson)+0.05*max(f_s_simpson),label=simletters[3])) +
  geom_text(aes(x=4, y=max(o_s_simpson)+0.05*max(o_s_simpson),label=simletters[4])) 
  pino2<-ggplot(shadata, aes(x=Substrate, y=Diversity)) +
  xlab("Sample")+ylab("Shannon Diversity")+ 
  geom_boxplot(outlier.shape=NA,fill=mycol) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0)) +
  geom_text(aes(x=1, y=max(f_p_shannon)+0.05*max(f_p_shannon),label=shaletters[1])) +
  geom_text(aes(x=2, y=max(o_p_shannon)+0.05*max(o_p_shannon),label=shaletters[2])) +
  geom_text(aes(x=3, y=max(f_s_shannon)+0.05*max(f_s_shannon),label=shaletters[3])) +
  geom_text(aes(x=4, y=max(o_s_shannon)+0.05*max(o_s_shannon),label=shaletters[4])) 
grid.arrange(pino,pino2,nrow=1)
dev.off()


shannon <- diversity(mydata)


#
limit<-min(rowSums(mydata))
simpson <- diversity(mydata, "simpson")
shannon <- diversity(mydata)
bray = vegdist(mydata, "bray") 
gower = vegdist(mydata, "gower")
pdf(paste(outdir,"rarefaction.pdf",sep=""))
rarecurve(mydata, col = "blue")
dev.off()
browser()
}
plot.diversity(infile=infile,outdir=outdir,rank=rank)


