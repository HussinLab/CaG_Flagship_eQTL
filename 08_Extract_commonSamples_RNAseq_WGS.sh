#Extract from WGS the samples in common with the RNAseq samples

mkdir WGS_common
cd WGS_common

#Get the list of common samples
#Put sample names of the WGS samples in the file WGS_samples.txt
grep -Ff <(cut -f 2 <path_to_metadata>/RNAseq.files.runs.911.noBadSamples.txt | sed '1d') WGS_samples.txt > WGS_ids_common_RNAseqQCed.txt
awk '{print $1"\t"$1}' WGS_ids_common_RNAseqQCed.txt > WGS_ids_common_RNAseqQCed.2cols.txt


#MAF 1% no LD pruning applied

for f in $(seq 1 22) ; do echo '#!/bin/bash' > chr$f.extract_and_filters.504Samples.maf0.01.noLDpruning.sh
echo "module load StdEnv/2020 plink/1.9b_6.21-x86_64
plink --vcf /lustre06/project/6065672/shared/Cartagene/flagship_paper/WGS/vcfs_freeze1_20221206/chr${f}_FINAL.vcf.gz --double-id --make-bed --real-ref-alleles --keep WGS_ids_common_RNAseqQCed.2cols.txt --out chr$f.matchRNAseq.PLINK.504
plink --bfile chr$f.matchRNAseq.PLINK.504 --geno 0.01 --make-bed --out chr$f.matchRNAseq.PLINK.504.miss1Perc --real-ref-alleles 
plink --bfile chr$f.matchRNAseq.PLINK.504.miss1Perc --maf 0.01 --make-bed --out chr$f.matchRNAseq.PLINK.504.miss1Perc.maf0.01 --real-ref-alleles 
plink --bfile chr$f.matchRNAseq.PLINK.504.miss1Perc.maf0.01 --hwe 0.000001 --make-bed --out chr$f.matchRNAseq.PLINK.504.miss1Perc.maf0.01.hwe0.000001 --real-ref-alleles " >> chr$f.extract_and_filters.504Samples.maf0.01.noLDpruning.sh

sbatch --chdir $(pwd) --account=ctb-hussinju --time=1:00:00 --output=chr$f.extract_and_filters.504Samples.maf0.01.noLDpruning.%j.out --error=chr$f.extract_and_filters.504Samples.maf0.01.noLDpruning.%j.err --mem=10G --cpus-per-task=12 chr$f.extract_and_filters.504Samples.maf0.01.noLDpruning.sh ; done

#merge autosomal
mkdir merged_maf1_noLDpruning
for f in $(seq 1 22) ; do echo "chr$f.matchRNAseq.PLINK.504.miss1Perc.maf0.01.hwe0.000001" >> merged_maf1_noLDpruning/files.txt ; done

module load StdEnv/2020
module load plink/1.9b_6.21-x86_64
plink --merge-list merged_maf1_noLDpruning/files.txt --make-bed --out merged_maf1_noLDpruning/autosomal.miss1Perc.maf0.01.hwe0.000001 --real-ref-alleles

#Remove low complexity from HG38.Anderson2010.bed definition
awk '{print $0"\t"NR}' HG38.Anderson2010.bed  | sed 's/^chr//' > lowComplexity.2.noCHR.bed

plink --bfile autosomal.miss1Perc.maf0.01.hwe0.000001 --exclude range lowComplexity.2.noCHR.bed--make-bed --real-ref-alleles --out autosomal.miss1Perc.maf0.01.hwe0.000001.woLowComplexity
