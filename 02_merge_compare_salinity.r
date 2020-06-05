library(data.table)
library(metagenomeSeq)

write.test.abundances<-function(minreads=1000,analysis="Peressotti")
{
#analysis<-"Peressotti"
#analysis<-"Basic_35_vs_0"

merged_seqs<-fread("basepath/tables/16s_merged_nochim.txt",data.table=F)
merged_taxa<-fread("basepath/tables/Taxa_merged_species.txt",data.table=F)
full_merged<-merge(merged_seqs,merged_taxa,by="V1",sort=F)
perc<-FALSE
if(perc)
{
for(aaa in 1:17)
{
full_merged[,as.character(aaa)]<-round(100*full_merged[,as.character(aaa)]/sum(full_merged[,as.character(aaa)]),5)
}
}


#Set up a metadata on the fly
if(analysis=="Basic_35_vs_0")
{
Group<-c(35,35,35,18,18,18,9,9,9,4,4,4,0,0,0,"FTQ","FD")
sample<-seq(1,17)
metadata<-data.frame(sample,Group,check.names=F)
#Set up the comparison (here 35psu vs 0)
tocompare<-c(35,0)
samplecomp<-metadata$sample[metadata$Group%in%tocompare]
metadata<-metadata[metadata$Group%in%tocompare,]
}

#We chose this one (suggested by Peressotti), which we think more interesting, contrasting 18 and 35 vs all the lower values of salinity.
#This is also supported by clustering results
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



full_merged$Species<-paste(full_merged$Genus,full_merged$Species,sep="_")

for (tax in c("Phylum","Class","Order","Family","Genus","Species"))
{
bytaxa<-aggregate(full_merged[,as.character(seq(1,17))],by=list(full_merged[,tax]),FUN="sum")
row.names(bytaxa)<-bytaxa$Group.1
bytaxa$Group.1<-NULL
bytaxa$Tot<-rowSums(bytaxa)
bytaxa<-bytaxa[order(bytaxa$Tot,decreasing=T),]
if(perc) mytype<-"perc" else mytype<-"count"

write.table(bytaxa,paste("basepath/tables/16s_merged_",mytype,"_by_",tax,".txt",sep=""),quote=F,col.names=NA,sep="\t")
#if(tax=="Order") browser()


#The analysis below make only sense if abundances are normalized i.e. expressed as percentages or reads per million 
if(!perc)
{
#Samples 16 and 17 are "special" samples and should be discarded
bytaxa$"16"<-bytaxa$"17"<-NULL
bytaxa<-bytaxa[bytaxa$Tot>=minreads,]
bytaxa$Tot<-NULL
mycomp<-bytaxa
#Only select samples that we want to compare
mycomp<-mycomp[,as.character(samplecomp)]
metad<-metadata
row.names(metad)<-metad$sample
metad$sample<-NULL
#Several steps below are needed to run the fitZig
#Create MRexperiment object
myciccio<-AnnotatedDataFrame(metad)
myexp<-newMRexperiment(counts=mycomp,phenoData=myciccio)
#Normalize
p <- cumNormStatFast(myexp)
Nexp <- cumNorm(myexp, p = p)
# Export count matrix
#Nexp <- MRcounts(Nexp, norm = TRUE, log = TRUE)
status <- as.character(pData(Nexp)$Group)
status <- factor(status)
mod <- model.matrix(~ status)
settings <- zigControl(maxit = 10, verbose = TRUE)
fit <- fitZig(obj = Nexp, mod = mod, useCSSoffset = FALSE, control = settings)
coefs <- MRcoefs(fit, coef = 2, group = 3, number = 1000)
final<-merge(coefs,mycomp,by="row.names")
final<-final[order(final$adjPvalues),]
#Row names of coefs are the sequences IDs.
write.table(final,paste("basepath/tables/test_differences/",analysis,"_by_",tax,".txt",sep=""),quote=F,sep="\t",row.names=F)
}
}

}

















