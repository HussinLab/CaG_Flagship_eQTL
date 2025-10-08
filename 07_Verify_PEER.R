library(ggplot2)
residuals=read.table("residuals_peer_60_allGenes.21082024.woOutliers.txt",h=T)
pca=prcomp(t(residuals))

pca_x=pca$x
pca_x=as.data.frame(pca_x)

mat=read.table("<path_to_data>/Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.txt",h=T,row.names=1)

toRem=c("11130732","11102119")

batchInfos=read.table("<path_to_metadata>/RNAseq.files.runs.911.txt",h=T)
batchInfos_2=as.character(batchInfos$Batch)
names(batchInfos_2)=batchInfos$Sample
names(mat)=gsub("X","",names(mat))
names(residuals)=names(mat)[!names(mat) %in% toRem]

batchInfos_911=as.factor(batchInfos_2[as.character(names(residuals))])
batchInfos_911_df=as.data.frame(batchInfos_911)
pca_x$Batch=as.factor(batchInfos_2[as.character(names(residuals))])
row.names(pca_x) = names(residuals)

pdf("pca.rnaseq.cag.allGenes.Combat.normalized.newFilters.residualsPeer60.21082024.noOutliers.pdf")
ggplot(pca_x,aes(PC1,PC2,col=Batch))+geom_point()
ggplot(pca_x,aes(PC3,PC4,col=Batch))+geom_point()
ggplot(pca_x,aes(PC5,PC6,col=Batch))+geom_point()
ggplot(pca_x,aes(PC7,PC8,col=Batch))+geom_point()
ggplot(pca_x,aes(PC9,PC10,col=Batch))+geom_point()
dev.off()

