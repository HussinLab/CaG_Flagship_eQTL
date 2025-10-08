#Look for any outlier after PEER first correction

library(ggplot2)
residuals=read.table("residuals_peer_60_allGenes.21082024.txt",h=T)

#Do a PCA with PEERs residuals
pca=prcomp(t(residuals))
pca_x=pca$x
pca_x=as.data.frame(pca_x)

mat=read.table("<path_to_normLCPM_counts>/Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.txt",h=T,row.names=1)

batchInfos=read.table("<path_to_metadata>/RNAseq.files.runs.911.txt",h=T)
batchInfos_2=as.character(batchInfos$Batch)
names(batchInfos_2)=batchInfos$Sample

#Add informations to the PEER table
names(residuals)=names(mat)
names(residuals)=gsub("X","",names(residuals))

#Add more informations about the batches
batchInfos_911=as.factor(batchInfos_2[as.character(names(residuals))])
batchInfos_911_df=as.data.frame(batchInfos_911)
pca_x$Batch=as.factor(batchInfos_2[as.character(names(residuals))])
row.names(pca_x) = names(residuals)


#Visualize the PCA
pdf("pca.rnaseq.cag.allGenes.Combat.normalized.newFilters.residualsPeer60.21082024.pdf")
ggplot(pca_x,aes(PC1,PC2,col=Batch))+geom_point()
ggplot(pca_x,aes(PC3,PC4,col=Batch))+geom_point()
ggplot(pca_x,aes(PC5,PC6,col=Batch))+geom_point()
ggplot(pca_x,aes(PC7,PC8,col=Batch))+geom_point()
ggplot(pca_x,aes(PC9,PC10,col=Batch))+geom_point()
dev.off()




#Look at the top loadings for PC1 and PC2
pca_loadings=as.data.frame(pca$rotation)
row.names(pca_loadings)=row.names(mat)

head(pca_loadings[order(pca_loadings$PC1,decreasing=T),c("PC1","PC2")])
#ENSG00000204520
head(pca_loadings[order(pca_loadings$PC2,decreasing=T),c("PC1","PC2")])
#ENSG00000204287

#FOUND 2 outliers 
> subset(pca_x, PC1 > 10)[,c(1,2,912)]
#              PC1       PC2  Batch
#11130732 56.43093 -1.491329 run481
> subset(pca_x, PC2 > 10)[,c(1,2,912)]
#              PC1      PC2  Batch
#11102119 3.464312 19.13366 run510


###### END OF PEER ROUND 1 ######