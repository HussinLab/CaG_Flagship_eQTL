####PEER CORRECTION :

library(peer)

#Load table
expression=read.table("/data/Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.txt",h=T,row.names="GeneID")
expression=expression[,1:length(names(expression))]
names(expression)=gsub("X","",names(expression))


covs = read.table('/data/RNAseq.files.runs.911.txt',header=TRUE)
covs_2=as.character(covs$Batch)
names(covs_2)=covs$Sample

covs_480=as.factor(covs_2[as.character(names(expression))])

covs_480_df=as.data.frame(covs_480)

#Initialize the model
model=PEER()

PEER_setPhenoMean(model,as.matrix(expression))
#because of our sample size : 60
PEER_setNk(model,60)
PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
pdf("/data/peer_correction/Precision_peer_60_allGenes.21082024.pdf")
plot(precision)
dev.off()

write.table(residuals,file="/data/peer_correction/residuals_peer_60_allGenes.21082024.txt")
write.table(factors,file="/data/peer_correction/factors_peer_60_allGenes.21082024.txt")
write.table(weights,file="/data/peer_correction/weights_peer_60_allGenes.21082024.txt")
write.table(precision,file="/data/peer_correction/precision_peer_60_allGenes.21082024.txt")

