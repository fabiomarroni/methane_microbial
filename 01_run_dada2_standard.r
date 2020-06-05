#This version is based on the standard dada2 workflow (https://benjjneb.github.io/dada2/tutorial.html)


library(dada2)

path <- "basepath/raw_reads"
ncores=8

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: XXX-YYY-SAMPLENAME_ZZZ.fastq
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 4)
sample.names <- sapply(strsplit(sample.names, "_"), `[`, 1)

#Plot some quality profile 

pdf("basepath/plots/quality_F.pdf")
plotQualityProfile(fnFs[c(1,13)])
dev.off()

pdf("basepath/plots/quality_R.pdf")
plotQualityProfile(fnRs[c(5,16)])
dev.off()

filtFs <- gsub("raw_reads","trimmed_reads", fnFs)
filtRs <- gsub("raw_reads","trimmed_reads", fnRs)
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#1) PRimers used for 16s amplification are:
# CCTACGGGNBGCASCAG
# GACTACNVGGGTATCTAATCC
# Lenghts are 17 and 21, respectively
#We set truncLen to 0 so that we do not put thresholds on read lengths

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230),
              trimLeft=c(17,21),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=ncores) # On Windows set multithread=FALSE


head(out)

#TruncLen=0
#                                         reads.in reads.out
#ID1432-16S-A01-1_S1_L001_R1_001.fastq.gz   141076     80890
#ID1432-16S-A02-2_S2_L001_R1_001.fastq.gz   149491     85865
#ID1432-16S-A03-3_S3_L001_R1_001.fastq.gz   156049     91500
#ID1432-16S-A04-4_S4_L001_R1_001.fastq.gz   180692    106774
#ID1432-16S-A05-5_S5_L001_R1_001.fastq.gz   159613    100389
#ID1432-16S-A06-6_S6_L001_R1_001.fastq.gz   160829     94168

#TruncLen=c(280,270)
#                                         reads.in reads.out
#ID1432-16S-A01-1_S1_L001_R1_001.fastq.gz   141076    104660
#ID1432-16S-A02-2_S2_L001_R1_001.fastq.gz   149491    110553
#ID1432-16S-A03-3_S3_L001_R1_001.fastq.gz   156049    117294
#ID1432-16S-A04-4_S4_L001_R1_001.fastq.gz   180692    135515
#ID1432-16S-A05-5_S5_L001_R1_001.fastq.gz   159613    124443
#ID1432-16S-A06-6_S6_L001_R1_001.fastq.gz   160829    120595




#Compute new error probabilitites
errF <- learnErrors(filtFs, multithread=ncores)
errR <- learnErrors(filtRs, multithread=ncores)


pdf("basepath/plots/errors_F.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()
pdf("basepath/plots/errors_R.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#Prepare dada sets using the newly estimated errors on the filtered reads
dadaFs <- dada(filtFs, err=errF, multithread=ncores)
dadaRs <- dada(filtRs, err=errR, multithread=ncores)

#Merge Forward and reverse reads
#Merge Forward and reverse reads, but when they do not merge then concatenate them
#(see https://github.com/benjjneb/dada2/issues/279)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, returnRejects=TRUE)
concat <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=TRUE)
#After this I want to keep two different data frames, to check performance of the merged and not merged ASV.

#This will create different sets for merged and non-merged reads
newm<-mergers
newc<-concat
for (aaa in 1:length(dadaFs))
{
newm[[aaa]]<-mergers[[aaa]][mergers[[aaa]]$accept,]
newc[[aaa]]<-concat[[aaa]][!mergers[[aaa]]$accept,]
} 


#Build sequence tabs
merged_seqtab <- makeSequenceTable(newm)
concat_seqtab <- makeSequenceTable(newc)

#Remove chimeras
merged_seqtab.nochim <- removeBimeraDenovo(merged_seqtab, method="consensus", multithread=ncores, verbose=TRUE)
concat_seqtab.nochim <- removeBimeraDenovo(concat_seqtab, method="consensus", multithread=ncores, verbose=TRUE)
dim(concat_seqtab.nochim )
#Change the colnames of concat sequence to split the sequences were concatenation put the NNs.
#Adding species with concatenated reads is not possible. Thus we follow the suggestion of splitting the read at the Ns and only feed the first portion to the addSpecies tool

colnames(concat_seqtab.nochim)<-unlist(lapply(strsplit(colnames(concat_seqtab.nochim),"NNNNNNNNNN"),"[",1))


#This is a very important result, since it gives the count of each sequence per sample!
merged_nochim<-t(merged_seqtab.nochim)
write.table(merged_nochim,"basepath/tables/16s_merged_nochim.txt",quote=F,col.names = NA)
concat_nochim<-t(concat_seqtab.nochim)
write.table(concat_nochim,"basepath/tables/16s_concat_first_nochim.txt",quote=F,col.names = NA)

#If something go bad we try to reread
emergency_recover=FALSE
if(emergency_recover)
{
concat_nochim<-fread("basepath/tables/16s_concat_first_nochim.txt",data.table=F)
pino<-aggregate(concat_nochim[,2:ncol(concat_nochim)],by=list(concat_nochim$V1),FUN="sum")
concat_nochim<-pino
row.names(concat_nochim)<-concat_nochim$Group.1
concat_nochim$Group.1<-NULL
write.table(concat_nochim,"basepath/tables/16s_concat_first_unique_nochim.txt",quote=F,col.names = NA)
concat_seqtab.nochim<-t(concat_nochim)

merged_nochim<-fread("basepath/tables/16s_merged_nochim.txt",data.table=F)
pino<-aggregate(merged_nochim[,2:ncol(merged_nochim)],by=list(merged_nochim$V1),FUN="sum")
merged_nochim<-pino
row.names(merged_nochim)<-merged_nochim$Group.1
merged_nochim$Group.1<-NULL
merged_seqtab.nochim<-t(merged_nochim)


}


#Track how many reads passed the QC steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(merged_seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,"basepath/tables/16s_sequence_stats",quote=F,col.names = NA)


merged_taxa <- assignTaxonomy(merged_seqtab.nochim, "/projects/populus/ep/share/marroni/databases/silva/silva_nr_v132_train_set.fa", multithread=ncores)
write.table(merged_taxa,"basepath/tables/Taxa_merged.txt",quote=F,col.names = NA)
concat_taxa <- assignTaxonomy(concat_seqtab.nochim, "/projects/populus/ep/share/marroni/databases/silva/silva_nr_v132_train_set.fa", multithread=ncores)
write.table(concat_taxa,"basepath/tables/Taxa_concat.txt",quote=F,col.names = NA)
cat("Remember to try to assing species on th especies silva database!!\n")


merged_taxa.species <- addSpecies(merged_taxa, "/projects/populus/ep/share/marroni/databases/silva/silva_species_assignment_v132.fa")
#Adding species with concatenated reads is not possible. Thus we follow the suggestion of splitting the read at the Ns and only feed the first portion to the addSpecies tool
write.table(merged_taxa.species,"basepath/tables/Taxa_merged_species.txt",quote=F,col.names = NA)

concat_taxa.species <- addSpecies(concat_taxa, "/projects/populus/ep/share/marroni/databases/silva/silva_species_assignment_v132.fa")
write.table(concat_taxa.species,"basepath/tables/Taxa_concat_species.txt",quote=F,col.names = NA)

cat("In assenza di errori classificate tutte le reads di Gloria!\n")


