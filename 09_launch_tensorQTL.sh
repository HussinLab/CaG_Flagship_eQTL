#Launch tensorqtl



#MAKE TRANSCRIPTOMIC BED FILE SHOULD LOOK LIKE :  (first line currently starts with #)

#Chr    start   end     TargetID        11100008        11100029        11100053        11100081        11100136        11100147        11100154        11100342        11100449        11100509        11100811   >
#1       923928  923929  ENSG00000187634 2.05438664140421        2.26350918294004        1.9349719990632 2.47255090951024        2.13596893410527        2.12249094836623        2.07131241828142        1.796394101>
#1       959309  959310  ENSG00000188976 5.62545944012696        5.89287882267793        5.78035702919193        6.05071852866138        5.71952141876566        5.88697126993732        5.94197848394526        5.6>
#1       960584  960585  ENSG00000187961 3.81541761545498        3.2955514038905 3.00438380633538        3.32094983862895        3.84664115264984        3.85363196408817        3.08823417398146        3.358481497>
#1       1001138 1001139 ENSG00000187608 5.1375219026663 3.9714484224517 5.84850087367494


#TSS FILE : EXTRACT gene lines (3rd column == gene) and look if column 7 is + or -, then take tss accordingly (column 4 or 5)

####




cd <path_to_tensorQTL_working_dir>/input

cp <path_to_data>/Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.txt .
sed -e 's/Sample_//' -e 's/_MPS.*//' <path_to_metadata>/WGS_ids_common_RNAseqQCed.txt > matchWGS_RNAseqQCed.txt

cp <path_to_data>/autosomal.miss1Perc.maf0.01.hwe0.000001.{bed,bim,fam} .
cp autosomal.miss1Perc.maf0.01.hwe0.000001.fam autosomal.miss1Perc.maf0.01.hwe0.000001.fam.bkp

#change names in the raw files so it matches RNAseq ids
sed -i -e 's/Sample_//g' -e 's/_MPS\S\+//g' autosomal.miss1Perc.maf0.01.hwe0.000001.fam


#NEW R SESSION to reorder samples correctly in the RNAseq, to match with WGS samples and add genes informations needed to run tensorQTL to the INT Normalized file

echo "genes_t=as.data.frame(read.table(\"Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.txt\",h=T,row.names=1))
genes=t(genes_t)

#extract the good samples too
samples=read.table(\"matchWGS_RNAseqQCed.txt\",h=F)

tss=read.table(\"tss.gtf\",h=T)
#genes_tss=merge(genes,tss,by.x=\"GeneID\",by.y=\"gene\")
genes_tss=merge(genes,tss,by.x=\"row.names\",by.y=\"gene\")
genes_tss_504=genes_tss[,names(genes_tss) %in% c(\"Row.names\",samples\$V1,\"chr\",\"start\",\"end\")]

#reorder
sample_order=read.table(\"autosomal.miss1Perc.maf0.01.hwe0.000001.fam\",h=F)
genes_tss_504_rd=genes_tss_504[,c(1,match(sample_order\$V1,names(genes_tss_504)),506,507,508)]

write.table(genes_tss_504_rd[,c(506,507,508,1:505)],file=\"Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.bed\",quote=F,row.names=F)" > reorderRNAseq.R

Rscript reorderRNAseq.R

#Sort the file correctly with the genes coordinates genome wide
sort -k1,1n -k2,2n Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.bed > Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort.bed

#Need to put autosomal first and sexual chromosomes after
awk '{if($1=="MT"|| $1=="X" || $1=="Y")print}' Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort.bed | sort -k1,1 -k2,2n > tmp_xym.txt
cat <(awk '{if($1!="MT" && $1!="X" && $1!="Y")print}' Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort.bed) tmp_xym.txt > Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort2.bed

#rm temporary files
rm Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort.bed tmp_xym.txt Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.bed

#Need to finally remove alternative haplotypes, extract sexual and MT chromosomes to order them then merge them back... change spaces for tabs
sed 's/ /\t/g' Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort2.bed >Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort.bed
rm Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort2.bed

############################
#### !!!!!!!Manually MOD HEADER AND PUT TABS BEFORE ZIPPING !!!!!!#####
#Chr start end TargetID
############################

#GZIP
gzip Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort.bed

##########################
#CREATE the metadata table 
##########################


#Extract samples and reorder PEER :
#Remove outliers : 
cut -f1 ../Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.txt > samples.911.txt

#outliers will be in outliers.txt file 2 outliers
grep -Fwvf outliers.txt samples.911.txt > samples.noOutliers.909.txt
#this will be used to put sample names to the PEER values and added to the final metadata file

#ADD COLUMN IN HEADER of weights_peer_60_allGenes.21082024.woOutliers.txt FILE AND THEN REORDER
paste samples.noOutliers.909.txt <(sed 's/ /\t/g' <path_to_PEER_results>/weights_peer_60_allGenes.21082024.woOutliers.txt | cut -f 2-) > weights_peer_60_allGenes.21082024.woOutliers.wsampleNames.txt
sed -i '1s/V/PEER/g' weights_peer_60_allGenes.21082024.woOutliers.wsampleNames.txt



#R SESSION to reorder PEER in the same order as teh WGS samples
echo "peer_weights=read.table(\"weights_peer_60_allGenes.21082024.woOutliers.wsampleNames.txt\",h=T,row.names=1)
sample_order=read.table(\"autosomal.miss1Perc.maf0.01.hwe0.000001.fam\",h=F)
peer_weights_504_rd=t(peer_weights[match(sample_order\$V1,row.names(peer_weights)),])

write.table(peer_weights_504_rd,file=\"weights_peer_60_allGenes.21082024.504.tr.txt\",sep=\"\t\",quote=F)" > reorder_PEER_table.R 
Rscript reorder_PEER_table.R 




### Use sovariates from Taliun's lab : 
cp <path_to_GWAS_covariates>/GWAS_ALL_covariates.tsv .
head -n 1 GWAS_ALL_covariates.tsv > GWAS_ALL_covariates.504.tsv ; grep -Fwf <(cut -d ' ' -f 1 autosomal.miss1Perc.maf0.01.hwe0.000001.fam) GWAS_ALL_covariates.tsv >> GWAS_ALL_covariates.504.tsv

#Transpose the resulting covariates file
awk '{ for ( i=1; i <=NF;i++ ) row[i] = row[i]((row[i])?" ":"")$i}; END{ for ( x = 1; x <= length(row); x++ )print row[x]}' GWAS_ALL_covariates.504.tsv | sed -e '1d' -e 's/ /\t/g' > GWAS_ALL_covariates.504.tr.tsv
#change the header a bit
sed -i '1s/IID/id/' GWAS_ALL_covariates.504.tr.tsv



#Launch tensorqtl
#DO ANALYSIS ON MAF1% no LD

##launch with the same covariates as the GWAS and PHEWAS (10 PCs, genotyping array, Age, Sex, Sex^2), INT NORMALIZATION

cat GWAS_ALL_covariates.504.tr.tsv <(sed '1d' weights_peer_60_allGenes.21082024.504.tr.txt) > metadata.eQTLs.matchRNAseq_WGS.CaG.504.60PEER.10PCs.array.Age.Sex.Age2.txt

cd <path_to_tensorQTL_working_dir>
mkdir output_withHLA

#Create header
echo '#!/bin/bash' > header_tQTL.allGenes.INT.withHLA.sh
echo "module load StdEnv/2020  gcc/9.3.0  cuda/11.4
module load python/3.8.10
module load arrow/2.0.0
module load r/4.1.0
python -c \"import pyarrow\"

source <path_to_python_virtualenv>

cd <path_to_tensorQTL_working_dir>

input_genetics=input/autosomal.miss1Perc.maf0.01.hwe0.000001.woLowComplexity
input_transcriptomics=input/Matrix.CaG.GRCh38.E98.kallisto.genes.Combat.aveLCPM0.normLCPM.21082024.INT_Norm.withTSS.sort.bed.gz" >> header_tQTL.allGenes.INT.withHLA.sh 

#Create the full bash file to run tensorQTL
cat header_tQTL.allGenes.INT.withHLA.sh > tQTL.10PCs.Array.Age.Age2.Sex.PEER.INT.withHLA.sh
echo "python3 -m tensorqtl --mode cis --cis_output output_withHLA/cis_PEER_array_Age_sex_Age2.signifOnly_withHLA -o output_withHLA/ --covariates input/metadata.eQTLs.matchRNAseq_WGS.CaG.504.60PEER.10PCs.array.Age.Sex.Age2.txt \$input_genetics \$input_transcriptomics cis_signif_PEER_array_Age_sex_Age2_withHLA
python3 -m tensorqtl --mode cis_nominal -o output_withHLA/ --covariates input/metadata.eQTLs.matchRNAseq_WGS.CaG.504.60PEER.10PCs.array.Age.Sex.Age2.txt \$input_genetics \$input_transcriptomics cis_nominal_tensorQTL_PEER_array_Age_sex_Age2_withHLA
" >> tQTL.10PCs.Array.Age.Age2.Sex.PEER.INT.withHLA.sh

#USE GPUs
sbatch --chdir $(pwd) --account=ctb-hussinju --time=1:00:00 --gres=gpu:1 --output=tQTL.10PCs.Array.Age.Age2.Sex.PEER.INT.withHLA.%j.out --error=tQTL.10PCs.Array.Age.Age2.Sex.PEER.INT.withHLA.%j.err --cpus-per-task=28 --mem=178g tQTL.10PCs.Array.Age.Age2.Sex.PEER.INT.withHLA.sh


