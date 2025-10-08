#FOR GWAS COVARIATES and GWAS


###NOTES ON RNASEQ 504 composition (clusters from Simon's analysis)
    #   5 0.0 (Haitian)
    #  17 10.0 (French Canadian)
    #   2 11.0 (French Canadian)
    #  20 12.0 (French Canadian)
    # 114 13.0 (French Canadian)
    # 343 14.0 (French Canadian)
    #   3 4.0 (Moroccan)



cd /lustre06/project/6065672/shared/Cartagene/flagship_paper/eQTLs/work_on_911_SampleFilters_after/tensorQTL/output_hla/gwas_covariates



### WITH R


library(arrow)
library(qvalue)

chr1_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.1.parquet")
chr2_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.2.parquet")
chr3_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.3.parquet")
chr4_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.4.parquet")
chr5_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.5.parquet")
chr6_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.6.parquet")
chr7_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.7.parquet")
chr8_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.8.parquet")
chr9_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.9.parquet")
chr10_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.10.parquet")
chr11_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.11.parquet")
chr12_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.12.parquet")
chr13_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.13.parquet")
chr14_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.14.parquet")
chr15_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.15.parquet")
chr16_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.16.parquet")
chr17_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.17.parquet")
chr18_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.18.parquet")
chr19_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.19.parquet")
chr20_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.20.parquet")
chr21_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.21.parquet")
chr22_all_10pcs=read_parquet("../cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA.cis_qtl_pairs.22.parquet")

all_10pcs=rbind(chr1_all_10pcs,chr2_all_10pcs,chr3_all_10pcs,chr4_all_10pcs,chr5_all_10pcs,chr6_all_10pcs,chr7_all_10pcs,chr8_all_10pcs,chr9_all_10pcs,chr10_all_10pcs,chr11_all_10pcs,chr12_all_10pcs,chr13_all_10pcs,chr14_all_10pcs,chr15_all_10pcs,chr16_all_10pcs,chr17_all_10pcs,chr18_all_10pcs,chr19_all_10pcs,chr20_all_10pcs,chr21_all_10pcs,chr22_all_10pcs)




all_10pcs$variant_id=gsub("_m","",all_10pcs$variant_id)
all_10pcs$assoc=paste0(all_10pcs$phenotype_id,":",all_10pcs$variant_id)



all_10pcs_ids=as.data.frame(all_10pcs$variant_id)
names(all_10pcs_ids)="variant_id"

all_10pcs$fdr=qvalue(all_10pcs$pval_nominal)$qvalues


rm(chr1_all_10pcs)
rm(chr2_all_10pcs)
rm(chr3_all_10pcs)
rm(chr4_all_10pcs)
rm(chr5_all_10pcs)
rm(chr6_all_10pcs)
rm(chr7_all_10pcs)
rm(chr8_all_10pcs)
rm(chr9_all_10pcs)
rm(chr10_all_10pcs)
rm(chr11_all_10pcs)
rm(chr12_all_10pcs)
rm(chr13_all_10pcs)
rm(chr14_all_10pcs)
rm(chr15_all_10pcs)
rm(chr16_all_10pcs)
rm(chr17_all_10pcs)
rm(chr18_all_10pcs)
rm(chr19_all_10pcs)
rm(chr20_all_10pcs)
rm(chr21_all_10pcs)
rm(chr22_all_10pcs)




gtex_signif=read.table("/lustre06/project/6065672/shared/GTex/v8/cis_eQTLs/global/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",h=T)
gtex_signif$variant_id=gsub("_b38","",gtex_signif$variant_id)
gtex_signif$gene_id2=gsub("\\.\\d+","",gtex_signif$gene_id,perl=T)
gtex_signif$assoc=paste0(gtex_signif$gene_id2,":",gtex_signif$variant_id)


#12360 genes signif in GTEx out of 20059 genes tested

#8252 genes common with CaG


#use already generated file by tensorQTL to get the p-value threshold per gene for the FDR 5%

gunzip -c ../../cis_signif_PEER_array_Age_sex_Age2_withHLA.cis_qtl.txt.gz | cut -f 1,19 > pval_threshold_per_gene.txt


#FROM THE RSESSION ....
#OPEN AND LOOK FOR THAT FILE AND DEFINE FOR EACH GENE IF PASSED OR NOT

thresholds_per_gene=read.table("FDR_per_GENE/pval_threshold_per_gene.txt",h=T,row.names=1)

all_10pcs_thr_per_gene=merge(all_10pcs,thresholds_per_gene,by.x="phenotype_id",by.y="row.names")

all_10pcs_thr_per_gene$fdr5_per_gene="Negative"
all_10pcs_thr_per_gene[which(all_10pcs_thr_per_gene$assoc %in% gtex_signif$assoc),"fdr5_per_gene"]="Positive"
all_10pcs_thr_per_gene[which(all_10pcs_thr_per_gene$pval_nominal <= all_10pcs_thr_per_gene$pval_nominal_threshold & all_10pcs_thr_per_gene$fdr5_per_gene != "Positive"),"fdr5_per_gene"]="FDR5_only"
all_10pcs_thr_per_gene[which(all_10pcs_thr_per_gene$pval_nominal <= all_10pcs_thr_per_gene$pval_nominal_threshold & all_10pcs_thr_per_gene$fdr5_per_gene == "Positive"),"fdr5_per_gene"]="Positive_both"


all_10pcs_thr_per_gene$fdr5_per_gene=as.factor(all_10pcs_thr_per_gene$fdr5_per_gene)


write.table( all_10pcs_thr_per_gene[which(all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","FDR5_only")), ],file="FDR_per_GENE/all_chr.fdr5perc_perGene.tsv")

#summary(all_10pcs_thr_per_gene$fdr5_per_gene)
#FDR5_only      Negative      Positive Positive_both
#  2613269     100414934        368151        969174


#look for the merged...
#merge_gtex_CaG

merge_gtex_CaG_fdr_per_gene=merge(gtex_signif,all_10pcs_thr_per_gene,by="assoc")


write.table(merge_gtex_CaG_fdr_per_gene,file="FDR_per_GENE/Merged_Per_Assoc.GTEx_CaG.txt")


#LOOK FOR SOME STATS : 
#CARTAGENE  :
#Genes with a significant association
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","FDR5_only")),"phenotype_id"]))
# 12774

#eSNP significant
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","FDR5_only")),"variant_id"]))
# 2223590

#gene-variant significant associations :
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","FDR5_only")),"assoc"]))
# 4028514 



#GTEX  : 
#Genes with a significant association 
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","Positive")),"phenotype_id"]))
# 8346 (common genes signif gtex + tested) for which a gene has a same variant-gene tested in Cartagene and signif in gtex

#eSNP significant
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","Positive")),"variant_id"]))
# 968040

#gene-variant significant associations :
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","Positive")),"assoc"]))
# 1523724



#NUMBER OF GENES SIGNIF IN GTEX AND TESTED IN CAG : 
sum(unique(all_10pcs_thr_per_gene$phenotype_id) %in% unique(gtex_signif$gene_id2))
#[1] 8551


#NUMBER OF GENES SIGNIF IN CAG AND GTEX :
sum(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","FDR5_only")),"phenotype_id"]) %in% unique(gtex_signif$gene_id2))
#[1] 8452

#NUMBER OF GENES SIGNIF IN GTEX :
length(unique(gtex_signif$gene_id2))
#[1] 12360


#UNIQUE TO :

#CARTAGENE (only the assoc is unique, because FDR field was done for this one):
#Genes with a significant association
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("FDR5_only")),"phenotype_id"]))
# 12768

#eSNP significant
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("FDR5_only")),"variant_id"]))
# 1853567

#gene-variant significant associations :
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("FDR5_only")),"assoc"]))
# 2971066 


#GTEX  (only the assoc is unique, because FDR field was done for this one): 
#Genes with a significant association 
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive")),"phenotype_id"]))
# 7424

#eSNP significant
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive")),"variant_id"]))
# 374377

#gene-variant significant associations :
length(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive")),"assoc"]))
# 466276




write.table(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","FDR5_only")),"phenotype_id"]),file="FDR_per_GENE/Genes_signif_CaG.fdr5perc_perGene.txt",quote=F,row.names=F)
write.table(unique( all_10pcs_thr_per_gene[which( all_10pcs_thr_per_gene$fdr5_per_gene %in% c("Positive_both","Positive")),"phenotype_id"]),file="FDR_per_GENE/Genes_signif_GTEx.fdr5perc_perGene.txt",quote=F,row.names=F)



#in bash...

cd FDR_per_GENE

#WHERE TESTING IF AT LEAST ONE VARIANT-GENE IN COMMON WAS SIGNIFICANT IN GTEx (not necessarily in CaG)
grep -Fwf Genes_signif_GTEx.fdr5perc_perGene.txt Genes_signif_CaG.fdr5perc_perGene.txt > Genes_signif_commonCaG_GTEx.txt
# 8255 genes having same signif variant-gene in both

grep -Fwvf Genes_signif_commonCaG_GTEx.txt Genes_signif_GTEx.fdr5perc_perGene.txt > Genes_signif_GTEx_only.fdr5perc_perGene.txt
# 91 genes in GTEx only having same variant-gene tested

grep -Fwvf Genes_signif_commonCaG_GTEx.txt Genes_signif_CaG.fdr5perc_perGene.txt > Genes_signif_CaG_only.fdr5perc_perGene.txt
# 4519 genes in CaG only having same variant-gene tested





#WHERE THE GENE IS SIGNIFICANT BUT NOT FOR THE SAME VARIANT NECESSARILY

grep -Fwf ../../../output/gwas_covariates/FDR_per_GENE/Whole_Blood.v8.signif_eGenes.txt Genes_signif_CaG.fdr5perc_perGene.txt > Genes_signif_CaG_and_GTEx_notSameVariant.txt
# 8452 genes being significant in both


grep -Fwvf Genes_signif_CaG_and_GTEx_notSameVariant.txt ../../../output/gwas_covariates/FDR_per_GENE/Whole_Blood.v8.signif_eGenes.txt > Genes_signif_GTEx_only.fdr5perc_perGene.notSame_Variant_necessarily.txt
# 3909 genes

grep -Fwvf Genes_signif_CaG_and_GTEx_notSameVariant.txt Genes_signif_CaG.fdr5perc_perGene.txt > Genes_signif_CaG_only.fdr5perc_perGene.notSame_Variant_necessarily.txt
# 4323

# HOW MANY OF THOSE UNIQUE GENES TO CaG WERE TESTED IN GTEX
grep -Fwf Genes_signif_CaG_only.fdr5perc_perGene.notSame_Variant_necessarily.txt ../../../output/gwas_covariates/allTestedGenes_GTEx_Whole_Blood_V8.txt | wc -l
#3823





#LOOK WITH eQTL-GEN.....

#fixed that list, before, was taking all the tested genes, and not only the significant ones.
ln -s /lustre06/project/6065672/shared/eQTL-gen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.fdr5perc.ENSids.txt eQTL-gen_signifHits.ENSids.txt
ln -s /lustre06/project/6065672/shared/eQTL-gen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.ENSids.tested.txt eQTL-gen_tested.ENSids.txt


#WHEN CONSIDERING SAME GENE/VARIANT TESTED : 
grep -Fwvf eQTL-gen_signifHits.ENSids.txt Genes_signif_CaG_only.fdr5perc_perGene.txt > genes_noSignifGTEx_or_eQTL-gen.txt
#883 genes only signif in CaG and not in GTEx or eQTL-gen

grep -Fwf eQTL-gen_signifHits.ENSids.txt Genes_signif_CaG_only.fdr5perc_perGene.txt | wc -l
#3636 genes only signif in CaG and eQTL-gen, but not in GTEx


wc -l eQTL-gen_signifHits.ENSids.txt
#16987 eQTL-gen_signifHits.ENSids.txt signif in eQTL-gen...




#WHEN CONSIDERING A GENE SIGNIFICANT: 
grep -Fwvf eQTL-gen_signifHits.ENSids.txt Genes_signif_CaG_only.fdr5perc_perGene.notSame_Variant_necessarily.txt > genes_noSignifGTEx_or_eQTL-gen.notSameVariant_necessarily.txt
#847 genes only signif in CaG and not in GTEx or eQTL-gen

grep -Fwf eQTL-gen_signifHits.ENSids.txt Genes_signif_CaG_only.fdr5perc_perGene.notSame_Variant_necessarily.txt | wc -l
#3476 genes only signif in CaG and eQTL-gen, but not in GTEx

# HOW MANY ARE SIGNIFICANT IN CAG AND TESTED IN eQTL-GEN
grep -Fwf Genes_signif_CaG_only.fdr5perc_perGene.notSame_Variant_necessarily.txt eQTL-gen_tested.ENSids.txt | wc -l
#3885

# HOW MANY ARE SIGNIFICANT IN CAG AND TESTED IN GTEX, AND SIGNIFICANT IN EQTL-GEN
grep -Fwf <(grep -Fwf Genes_signif_CaG_only.fdr5perc_perGene.notSame_Variant_necessarily.txt ../../../output/gwas_covariates/allTestedGenes_GTEx_Whole_Blood_V8.txt) eQTL-gen_signifHits.ENSids.txt | wc -l
#3393

# HOW MANY ARE SIGNIFICANT IN CAG AND TESTED IN GTEX, AND TESTED IN EQTL-GEN
grep -Fwf <(grep -Fwf Genes_signif_CaG_only.fdr5perc_perGene.notSame_Variant_necessarily.txt ../../../output/gwas_covariates/allTestedGenes_GTEx_Whole_Blood_V8.txt) eQTL-gen_tested.ENSids.txt | wc -l
#3785

# HOW MANY ARE SIGNIFICANT IN GTEX ONLY AND TESTED IN CAG
grep -Fwf Genes_signif_GTEx_only.fdr5perc_perGene.notSame_Variant_necessarily.txt ../testedGenes.txt | wc -l
#99








cor.test(merge_gtex_CaG$slope_GTEx,merge_gtex_CaG$slope_CaG)

#        Pearson's product-moment correlation

#data:  merge_gtex_CaG$slope_GTEx and merge_gtex_CaG$slope_CaG
#t = 1780.4, df = 1337323, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.8381259 0.8391316
#sample estimates:
#      cor
#0.8386294
#Spearman :
#0.8514551




#best hits per gene per cohort

#GTEX :

library(dplyr)

gtex_signif_order=arrange(gtex_signif,gene_id2,pval_nominal)
gtex_signif_1hit=distinct(gtex_signif_order,gene_id2,.keep_all=1)

rm(gtex_signif_order)


#CaG
all_10pcs_thr_per_gene_order=arrange(all_10pcs_thr_per_gene,phenotype_id,pval_nominal)
all_10pcs_thr_per_gene_1hit=distinct(all_10pcs_thr_per_gene_order,phenotype_id,.keep_all=1)

rm(all_10pcs_thr_per_gene_order)


#MERGE BOTH
merged_1hit=merge(gtex_signif_1hit,all_10pcs_thr_per_gene_1hit,by.x="gene_id2",by.y="phenotype_id")


merged_1hit$diff=-log10(merged_1hit$pval_nominal.y) - -log10(merged_1hit$pval_nominal.x)


pdf("GTEx_vs_CaG.best_hit_perGene.slope_comparison.pdf")
#ggplot(merged_1hit,aes(x= slope.x, y=slope.y),color=diff)+geom_hex(bins=100)+xlim(c(-2.5,2.5))+ylim(c(-4,4))+xlab("Size Effect GTEx")+ylab("Size Effect Cartagene")
ggplot(merged_1hit,aes(x= slope.x, y=slope.y,color=diff,alpha=0.1))+geom_point()+xlim(c(-2.5,2.5))+ylim(c(-4,4))+xlab("Size Effect GTEx")+ylab("Size Effect Cartagene")+scale_color_gradient(low="yellow", high="purple4")
#ggplot(merged_1hit,aes(x= slope.x, y=slope.y),color=diff)+geom_hex(bins=50)+xlim(c(-1,1))+ylim(c(-2,2))+xlab("Size Effect GTEx")+ylab("Size Effect Cartagene")
dev.off()


pdf("GTEx_vs_CaG.best_hit_perGene.diff_pval.signif_both.pdf")
ggplot(subset(merged_1hit,pval_nominal.y <= pval_nominal_threshold.y) ,aes(x= diff))+geom_histogram(bins=50)+xlab("Difference between -log10 transformed Cartagene and GTEx pvalues")
dev.off()



dim(subset(merged_1hit,pval_nominal.y <= pval_nominal_threshold.y & pval_nominal.y < pval_nominal.x))
#[1] 5346   (out of 8452 common eGenes signif in CaG vs signif GTEx)

#8551 signif gtex vs not necessarily in cag

summary(subset(merged_1hit,pval_nominal.y < pval_nominal.x)$fdr5_per_gene)
#    FDR5_only      Negative      Positive Positive_both
#         2470             0             0          2876




write.table(subset(merged_1hit,pval_nominal.y < pval_nominal.x),file="FDR_per_GENE/Genes_better_bestHit_CaG_than_GTEx.txt")


write.table(merged_1hit,file="FDR_per_GENE/Better_hits_genes_common_GTEx_CaG.txt")










#GET THE EXTENT OF THE DIFFERENCE IN P-VALUES
#diff_best_in_cag=read.table("Genes_better_bestHit_CaG_than_GTEx.txt",h=T)
merged_1hit$diff=-log10(merged_1hit$pval_nominal.y) - -log10(merged_1hit$pval_nominal.x)



dim(subset(merged_1hit,pval_nominal.y <= pval_nominal_threshold.y & pval_nominal.y < pval_nominal.x))
#[1] 5346   29
dim(subset(merged_1hit,pval_nominal.y <= pval_nominal_threshold.y & pval_nominal.y < pval_nominal.x & diff < 10)) 
#[1] 2779   29 (or 2567 >= 10)
dim(subset(merged_1hit,pval_nominal.y <= pval_nominal_threshold.y & pval_nominal.y < pval_nominal.x & diff < 5))
#[1] 1705   29 (or 3641 >= 5)
dim(subset(merged_1hit,pval_nominal.y <= pval_nominal_threshold.y & pval_nominal.y < pval_nominal.x & diff < 2))
#[1] 821  29 (or 4525 >= 2)
dim(subset(merged_1hit,pval_nominal.y <= pval_nominal_threshold.y & pval_nominal.y < pval_nominal.x & diff < 1))
#[1] 436  29 (or 4910 >= 1 )









#COMPARE WITH ALL THE GTEX ASSOCIATIONS:
library(data.table)
gtex_all=as.data.frame(fread("/lustre06/project/6065672/shared/GTex/v8/cis_eQTLs/b37/Whole_Blood.allpairs.withB38_ids.txt"))
#gtex_all$variant_id=gsub("_b38","",gtex_all$variant_id)
gtex_all$gene_id2=gsub("\\.\\d+","",gtex_all$gene_id,perl=T)
gtex_all$assoc=paste0(gtex_all$gene_id2,":",gtex_all$variant_id_b38)


merge_gtex_all_CaG_fdr_per_gene=merge(gtex_all,all_10pcs_thr_per_gene,by="assoc")
dim(merge_gtex_all_CaG_fdr_per_gene)
#[1] 70347217       20
length(unique(merge_gtex_all_CaG_fdr_per_gene$gene_id2))
#12455


names(merge_gtex_all_CaG_fdr_per_gene)[7]="slope_GTEx"
names(merge_gtex_all_CaG_fdr_per_gene)[17]="slope_CaG"

cor.test(merge_gtex_all_CaG_fdr_per_gene$slope_GTEx,merge_gtex_all_CaG_fdr_per_gene$slope_CaG)


#         Pearson's product-moment correlation

# data:  merge_gtex_all_CaG_fdr_per_gene$slope_GTEx and merge_gtex_all_CaG_fdr_per_gene$slope_CaG
# t = 2468.8, df = 70347215, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2821581 0.2825882
# sample estimates:
#       cor
# 0.2823732

cor.test(merge_gtex_all_CaG_fdr_per_gene$slope_GTEx,merge_gtex_all_CaG_fdr_per_gene$slope_CaG,method="spearman")
#         Spearman's rank correlation rho

# data:  merge_gtex_all_CaG_fdr_per_gene$slope_GTEx and merge_gtex_all_CaG_fdr_per_gene$slope_CaG
# S = 4.9388e+22, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.1488078

# Warning message:
# In cor.test.default(merge_gtex_all_CaG_fdr_per_gene$slope_GTEx,  :
#   Cannot compute exact p-value with ties



pdf("GTEx_vs_CaG.slope_comparison.all_hits_GTEx.pdf")
ggplot(merge_gtex_all_CaG_fdr_per_gene,aes(x= slope_GTEx, y=slope_CaG))+geom_hex(bins=100)+xlim(c(-2.5,2.5))+ylim(c(-4,4))
ggplot(merge_gtex_all_CaG_fdr_per_gene,aes(x= slope_GTEx, y=slope_CaG))+geom_hex(bins=50)+xlim(c(-1,1))+ylim(c(-2,2))
dev.off()



#only signif GTEx :

gtex_signif=read.table("/lustre06/project/6065672/shared/GTex/v8/cis_eQTLs/global/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz",h=T)
gtex_signif$variant_id=gsub("_b38","",gtex_signif$variant_id)
gtex_signif$gene_id2=gsub("\\.\\d+","",gtex_signif$gene_id,perl=T)
gtex_signif$assoc=paste0(gtex_signif$gene_id2,":",gtex_signif$variant_id)


merge_gtex_CaG=merge(gtex_signif,all_10pcs,by="assoc")
dim(merge_gtex_CaG)
#[1] 1523724      23

library(tidyr)
names(merge_gtex_CaG)[8]="pval_nominal_GTEx"
names(merge_gtex_CaG)[21]="pval_nominal_CaG"


names(merge_gtex_CaG)[9]="slope_GTEx"
names(merge_gtex_CaG)[22]="slope_CaG"


pdf("GTEx_vs_CaG.slope_comparison.signifGTEx.pdf")
ggplot(merge_gtex_CaG,aes(x= slope_GTEx, y=slope_CaG))+geom_hex(bins=100)+xlim(c(-2.5,2.5))+ylim(c(-4,4))
ggplot(merge_gtex_CaG,aes(x= slope_GTEx, y=slope_CaG))+geom_hex(bins=50)+xlim(c(-1,1))+ylim(c(-2,2))
dev.off()



cor.test(merge_gtex_CaG$slope_GTEx,merge_gtex_CaG$slope_CaG,method="pearson")

#         Pearson's product-moment correlation

# data:  merge_gtex_CaG$slope_GTEx and merge_gtex_CaG$slope_CaG
# t = 1702, df = 1523722, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8089557 0.8100504
# sample estimates:
#       cor
# 0.8095038

cor.test(merge_gtex_CaG$slope_GTEx,merge_gtex_CaG$slope_CaG,method="spearman")


#         Spearman's rank correlation rho

# data:  merge_gtex_CaG$slope_GTEx and merge_gtex_CaG$slope_CaG
# S = 1.0636e+17, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.8196062

# Warning message:
# In cor.test.default(merge_gtex_CaG$slope_GTEx, merge_gtex_CaG$slope_CaG,  :
#   Cannot compute exact p-value with ties