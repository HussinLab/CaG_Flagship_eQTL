#START R SESSION
toRem=c("11130732","11102119")
expression=read.table("/data/Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.txt",h=T,row.names="GeneID")
expression=expression[,1:length(names(expression))]
names(expression)=gsub("X","",names(expression))
library(peer)

expression=expression[,!names(expression) %in% toRem]
model=PEER()

#PEER_setCovariates(model, as.matrix(covs_480_df))


PEER_setPhenoMean(model,as.matrix(expression))
#PEER_setNk(model,10)
PEER_setNk(model,60)
PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)
pdf("/data/peer_correction/Precision_peer_60_allGenes.21082024.woOutliers.pdf")
plot(precision)
dev.off()


write.table(residuals,file="/data/peer_correction/residuals_peer_60_allGenes.21082024.woOutliers.txt")
write.table(factors,file="/data/peer_correction/factors_peer_60_allGenes.21082024.woOutliers.txt")
write.table(weights,file="/data/peer_correction/weights_peer_60_allGenes.21082024.woOutliers.txt")
write.table(precision,file="/data/peer_correction/precision_peer_60_allGenes.21082024.woOutliers.txt")


