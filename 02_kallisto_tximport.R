library(GenomicFeatures)
library(tximport)
library(rhdf5)
tx2db=makeTxDbFromGFF("<path_to_gtf>/Homo_sapiens.GRCh38.98.gtf",format="gtf")
dir <- system.file("extdata", package = "tximportData")
k <- keys(tx2db, keytype = "TXNAME")
tx2gene <- select(tx2db, k, "GENEID", "TXNAME")

#need to put the name of the output folder of kallisto runs in the same file
samples <- read.table("<path_to_kallisto_results_rssource_file>/files.txt", header = TRUE)
files <- file.path(samples$Sample, "abundance.h5")
names(files)=gsub(".merged","",gsub(".GRCh38.*","",gsub(".*recoded.","",files,perl=T),perl=T))
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

txi.sum <- summarizeToGene(txi.kallisto, tx2gene,ignoreTxVersion=T)
write.table(txi.sum$counts,file="Matrix.911samples.RNAseq.kallisto.genes.unstranded.moreTranscripts.txt",quote=F)

txi.sum.lengthscaledtpm <- summarizeToGene(txi.kallisto, tx2gene_final,countsFromAbundance="lengthScaledTPM",ignoreTxVersion=T)
#transcripts missing from tx2gene: 18143

write.table(txi.sum.lengthscaledtpm$abundance,file="Matrix.911samples.RNAseq.kallisto.genes.unstranded.TPM.txt",quote=F)
