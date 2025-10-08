#1) Work on the RNAseq data, do all the QCs to know the remaining set of samples in common after with the WGS

#with R
library(edgeR)
library(sva)
library(ggplot2)

#Convert a GTF file to get a bed file for each genes first and import:
genes=read.table("<path_to_bed>/Homo_sapiens.GRCh38.109.genes.bed",h=F)

#import table
mat=read.table("<path_to_table>/Matrix.911samples.RNAseq.kallisto.genes.unstranded.txt",h=T,row.names=1)
names(mat)=gsub("^X","",names(mat),perl=T)


#Adjust for known batches with Combat
batchInfos=read.table("<path_to_metadata>/RNAseq.files.runs.911.txt",h=T)
batchInfos_2=as.character(batchInfos$Batch)
names(batchInfos_2)=batchInfos$Sample

dge=DGEList(mat,group=NULL)
mat_batch_adjusted <- ComBat_seq(as.matrix(mat), batch=as.character(batchInfos$Batch), group=NULL)
mat_batch_adjusted=as.data.frame(as.matrix(mat_batch_adjusted))

#New DGE object
dge=DGEList(mat_batch_adjusted,group=NULL)


#Filter lowly expressed genes using raw counts
average_lcpm <- aveLogCPM(dge)
summary(average_lcpm > 0)
#   Mode   FALSE    TRUE
#logical   21987   13641

completeGeneList=row.names(dge)
passingLCPM=completeGeneList[average_lcpm > 0]

dge_kept <- dge[passingLCPM , , keep.lib.sizes=FALSE]
#[1] 13641   

#Show the effect of normalization
raw_dge <- dge_kept
norm_dge <- calcNormFactors(dge_kept, method = "TMM")
raw_lcpm <- cpm(raw_dge, log = TRUE)
norm_lcpm <- cpm(norm_dge, log = TRUE)

norm_lcpm_withinfos=merge(norm_lcpm,genes,by.x="row.names",by.y="V4")

#Redefine the order of the columns
norm_lcpm_withinfos=norm_lcpm_withinfos[,c(913,914,915,1,2:912)]

#Get metadata from CaG and link them to the 111 ids
correspondance=read.table("<path_to_metadata>/correspondance.samples.RNAseq.<project_id>_ids.911.21082024.txt",h=T)
sex_info=read.csv("<path_to_metadata>/111_data_baseline.mod.6cols.csv",h=T)

#Extract the 911 samples and merge the samples informations to get information about sex in the metadata
sex_info_merge=merge(sex_info,correspondance,by="project_code")
sex_info_2=as.character(sex_info_merge$SEXE)
names(sex_info_2)=sex_info_merge$files_111
sex_infos_911=as.factor(sex_info_2[as.character(names(norm_lcpm_withinfos[5:915]))])
sex_infos_911_df=as.data.frame(sex_infos_911)


#Do a PCA and see the effects due to sex or if anything left for the Batches
pca=prcomp(t(norm_lcpm_withinfos[,5:915]))
pca_x=pca$x
pca_x=as.data.frame(pca_x)
pca_x$Batch=as.factor(batchInfos_2[names(norm_lcpm_withinfos[5:915])])
pca_x$Sex=as.factor(sex_infos_911[as.character(names(norm_lcpm_withinfos[5:915]))])

pdf("<path_to_results>/pca.rnaseq.cag.allGenes.Combat.normalizedLCPM.newFilters_LCPM.21082024.911.pdf")
ggplot(pca_x,aes(PC1,PC2,col=Batch))+geom_point()
ggplot(pca_x,aes(PC1,PC2,col=Sex))+geom_point()
dev.off()


#Sex mismatches identification :
#3 samples seems to don't have the good sex label (manual curation):
#Males are on the negative side, Female on the positive site of PC2
subset(pca_x, (PC2 > 10 & Sex == "MALE") | (PC2 < -5 & Sex == "FEMALE"))[,c("PC2","Batch","Sex")]
#                PC2          Batch    Sex
#11129736 -18.033310 run147-148-159 FEMALE
#11121123  13.695892 run147-148-159   MALE
#11118427  -6.522931         run589 FEMALE

#From them, 11121123 present in wgs 
samples_to_remove=c("11129736","11121123","11118427")

#Plot PCA with residuals coming from correction for the sex 
new.x = apply(as.matrix(norm_lcpm_withinfos[,5:915]), 1, FUN = function(x){return(resid(lm(as.numeric(x) ~ +as.factor(sex_infos_911[as.character(names(norm_lcpm_withinfos[5:915]))]) , na.action=na.exclude)))})
rownames(new.x)=colnames(norm_lcpm_withinfos[,5:915])


pca=prcomp(new.x)
pca_x=pca$x
pca_x=as.data.frame(pca_x)
pca_x$Batch=as.factor(batchInfos_2[as.character(names(norm_lcpm_withinfos[5:915]))])
pca_x$Sex=as.factor(sex_infos_911[as.character(names(norm_lcpm_withinfos[5:915]))])

pdf("pca.rnaseq.cag.allGenes.Combat.normalizedLCPM.newFilters_LCPM.SexCorrected.21082024.pdf")
ggplot(pca_x,aes(PC1,PC2,col=Batch))+geom_point()
ggplot(pca_x,aes(PC1,PC2,col=Sex))+geom_point()
dev.off()



#WRITE TABLE WITH BATCH NORMALIZATION AND LCPM VALUES
write.table(norm_lcpm,file="Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.txt",sep="\t")




#Normalization steps :
#Quantile Normalization :
qn <- function(.data){
 data_sort <- apply(.data, 2, sort)
 row_means <- rowMeans(data_sort)
 data_sort <- matrix(row_means, 
                     nrow = nrow(data_sort), 
                     ncol = ncol(data_sort), 
                     byrow = TRUE
                     )
 index_rank <- apply(.data, 2, order)
 normalized_data <- matrix(nrow = nrow(.data), ncol = ncol(.data))
 for(i in 1:ncol(.data)){
   normalized_data[,i] <- data_sort[index_rank[,i], i]
 }
 return(normalized_data)
}
normalized_data <- qn(norm_lcpm)
colnames(normalized_data)=colnames(norm_lcpm)
rownames(normalized_data)=rownames(norm_lcpm)

#WRITE TABLE WITH BATCH NORMALIZATION AND LCPM VALUES QUANTILE NORMALIZED
write.table(normalized_data,file="Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.quantileNorm.txt",sep="\t")

### MANUALLY ADD "GeneID" in the header of the created file


#INT Normalization :
library("RNOmni")
int_data <- apply(norm_lcpm,1,RankNorm)

#WRITE TABLE WITH BATCH NORMALIZATION AND LCPM VALUES INT NORMALIZED
write.table(int_data,file="Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.txt",sep="\t")
