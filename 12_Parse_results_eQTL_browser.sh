# eQTL browser :
#git clone https://github.com/HussinLab/fivex_eQTL

#need to create a python virtual environment
#python -v venv virtualenv
#
#
#need to source after :
#source virtualenv/bin/activate
#need to install Jinja2==3.1.4

#was installed using conda instead of virtual environment
#needed to re-install flask and a compatible version of vue

#For development :
#need to run both ./run-development.sh file and npm run serve commands (front- and back- ends running at the same time)




#need to generate :
# /data
# gene.id.symbol.map.json.gz
# rsid.sqlite3.db
# 
# /data 
#



#For the input files :
#Generate the references:

#from the raw data on dbsnp #TODELETE:
#gunzip -c dbsnp151.GRCh38.tsv.gz | awk '{print $1"_"$2"_"$3"_"$4"\t"$5}' > dbsnp151.GRCh38.2cols.tsv



#Create a temporary folder in the data directory to create the files we need :
mkdir <path_to_eqtl_browser>/fivex/data/cartagene_ge 
cd <path_to_eqtl_browser>/fivex/data/cartagene_ge


#CaG list :
#Created from R scripts in which Cartagene dataset was all compiled. This is in reality only the concatenated table from tensorQTL with an additional column that can just be ignored.

cp <path_to_eqtl_results>/all_chr.fdr5perc_perGene.tsv .


#NEED THIS HEADER : 
#"phenotype_id" "variant_id" "start_distance" "af" "ma_samples" "ma_count" "pval_nominal" "slope" "slope_se" "assoc" "pval_nominal_threshold" "fdr5_per_gene"
#ADJUST WITH THE RIGHT COLUMNS

#depending on the format, we would want in the end a tabulated file without row numbers at the beginning
head -n1 all_chr.fdr5perc_perGene.tsv > header
sed '1d' all_chr.fdr5perc_perGene.tsv | cut -d ' ' -f2- > noheader
cat header noheader | sed -e 's/"//g' -e 's/ /\t/g' > all_chr.fdr5perc_perGene.tsv

rm header noheader

#extract positions from the list and keep the original order
cut -f 2 all_chr.fdr5perc_perGene.tsv | sed -e 's/^chr//' -e 's/_/\t/g'  > CaG.variants.inOrder.tsv
sed -i '1s/variant\tid/chr\tpos\tref\talt/' CaG.variants.inOrder.tsv

#define, in the same order if a variant is a SNP or INDEL
awk '{ type="SNP" ; if(length($3) > 1 || length($4) > 1 ) type="INDEL" ; print type}' CaG.variants.inOrder.tsv | sed -e '1s/INDEL/type/' > CaG.variants.inOrder.types.tsv




#GTEx list :

#From the significant GTEx results (v8 version) :
cp <path_to_GTEx_Analysis_v8_eQTL>/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz .

#remove _b38 so we have variants compatible with our version and remove the gene version
gunzip -c Whole_Blood.v8.signif_variant_gene_pairs.txt.gz | sed -e 's/_b38//' -e 's/\t\(ENS\S\+\)\.[0-9]\+\t/\t\1\t/' > GTEx_ge_blood.all.tsv.tmp
rm Whole_Blood.v8.signif_variant_gene_pairs.txt.gz

#extract positions and define the type after
cut -f 1 GTEx_ge_blood.all.tsv.tmp | sed -e 's/^chr//' -e 's/_/\t/g' > GTEx.blood.variants.inOrder.tsv
sed -i '1s/variant\tid/chr\tpos\tref\talt/' GTEx.blood.variants.inOrder.tsv

#define type
awk '{ type="SNP" ; if(length($3) > 1 || length($4) > 1 ) type="INDEL" ; print type}' GTEx.blood.variants.inOrder.tsv | sed -e '1s/INDEL/type/' > GTEx.blood.variants.inOrder.types.tsv





#CREATE RS ID DATABASE AND MERGE THE POSITIONS TO THE OTHER TABLES AFTER

module load StdEnv/2020
module load bcftools/1.11

head -n 1 CaG.variants.inOrder.tsv > CaG.allVariants.uniq.tsv



cat <(sed '1d' CaG.variants.inOrder.tsv) <(sed '1d' GTEx.blood.variants.inOrder.tsv) | sort -k1,1n -k2,2n | uniq >> CaG.allVariants.uniq.tsv
cat CaG.allVariants.uniq.tsv | awk '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' | sed '1d' > CaG.allVariants.uniq.bed

#extract the right information in the good order from dbSNP file
gunzip -c <path_to_dbsnp_data>/dbsnp151.GRCh38.00-All.vcf.gz | grep -vP "^#" | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$3}' > dbsnp151.GRCh38.tsv 

#index the file to do a query
bgzip dbsnp151.GRCh38.tsv
tabix -s 1 -b 2 -e 2 dbsnp151.GRCh38.tsv.gz

#do the query of CaG variants :
tabix -R CaG.allVariants.uniq.bed dbsnp151.GRCh38.tsv.gz > CaG.allVariants.uniq.dbSNP.annot.tsv


#do a backup, just in case
cp CaG.allVariants.uniq.dbSNP.annot.tsv CaG.allVariants.uniq.dbSNP.annot.tsv.bkp

#will be used to create the database
cat CaG.allVariants.uniq.dbSNP.annot.tsv  | awk '{split ($4,t,",");for(i in t){$4=t[i];print $0}}' | sed 's/ /\t/g' > CaG.allVariants.uniq.dbSNP.annot.onePerLine.tsv


#will be used later ...
awk '{print $1"_"$2"_"$3"_"$4"\t"$5}' CaG.allVariants.uniq.dbSNP.annot.onePerLine.tsv > CaG.allVariants.uniq.dbSNP.annot.2cols.tsv


#Do the database :

sqlite3 rsid.sqlite3.db

create table rsidTable(chrom text,pos int, ref text, alt text, rsid text);
.mode tabs
.import CaG.allVariants.uniq.dbSNP.annot.onePerLine.tsv rsidTable

CREATE INDEX idx_chrom_pos ON rsidTable (chrom, pos);
CREATE INDEX idx_rsid ON rsidTable (rsid);


#exit




#RECOVER RSIDs for CaG
sed 's/\t/_/g' CaG.variants.inOrder.tsv > CaG.variants.inOrder.1col.tsv

#RECOVER RSIDs for GTEx
sed 's/\t/_/g' GTEx.blood.variants.inOrder.tsv > GTEx.blood.variants.inOrder.1col.tsv




#Annotate, in the same order, the variants in CaG and GTEx with R

cag_variants=read.table("CaG.variants.inOrder.1col.tsv",h=T)
gtex_variants=read.table("GTEx.blood.variants.inOrder.1col.tsv",h=T)
annotations=read.table("CaG.allVariants.uniq.dbSNP.annot.2cols.tsv",h=F)
names(annotations)=c("chr_pos_ref_alt","id")

library(dplyr)

merged_cag=left_join(cag_variants,distinct(annotations,chr_pos_ref_alt,.keep_all=1),by="chr_pos_ref_alt")
merged_gtex=left_join(gtex_variants,distinct(annotations,chr_pos_ref_alt,.keep_all=1),by="chr_pos_ref_alt")


write.table(merged_cag,file="CaG.variants.inOrder.1col.annot.tsv",quote=F,sep="\t",row.names=F)
write.table(merged_gtex,file="GTEx.blood.variants.inOrder.1col.annot.tsv",quote=F,sep="\t",row.names=F)

#exit




#DO SOME OF THE FILES THAT WILL BE IMPORTED BY FIVEX


#needs to be transferred in : {DATA_DIR}/ebi_original/{STUDY}/ge/{STUDY}ge{TISSUE}.all.tsv.gz
#FORMAT NEEDED :
# 2: molecular_trait_id
# 3: chromosome
# 4: position (int)
# 5: ref
# 6: alt
# 7: variant (chr_pos_ref_alt)
# 8: ma_samples (int)
# 9: maf (float)
# 10: pvalue (float)
# 11: beta (float)
# 12: se (float)
# 13: type (SNP, INDEL, etc)
# 14: ac (allele count) (int)
# 15: an (total number of alleles = 2 * sample size) (int)
# 16: r2 (float)
# 17: molecular_trait_object_id
# 18: gene_id (ENSG#)
# 19: median_tpm (float)
# 20: rsid


#temporary unsorted files

#Cartagene :
paste <(cut -f 1 all_chr.fdr5perc_only.tsv | sed '1d') <(sed '1d' CaG.variants.inOrder.tsv) <(sed '1d' all_chr.fdr5perc_only.tsv | awk '{print $2"\t"$5"\t"$4"\t"$7"\t"$8"\t"$9}') <(sed '1d' CaG.variants.inOrder.types.tsv) <(sed '1d' all_chr.fdr5perc_only.tsv | awk '{an=($6/$4)+0.5 ; print $6"\t"int(an)"\t1\t"$1"\t"$1"\t0"}') <(cut -f 2 CaG.variants.inOrder.1col.annot.tsv | sed '1d') > Cartagene_ge_blood.all.unsorted.tsv

#GTEx :
paste <(cut -f 2 GTEx_ge_blood.all.tsv.tmp | sed '1d') <(sed '1d' GTEx.blood.variants.inOrder.tsv) <(sed '1d' GTEx_ge_blood.all.tsv.tmp | awk '{print $1"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9}') <(sed '1d' GTEx.blood.variants.inOrder.types.tsv) <(sed '1d' GTEx_ge_blood.all.tsv.tmp | awk '{an=($5/$6)+0.5 ; print $5"\t"int(an)"\t1\t"$2"\t"$2"\t0"}') <(cut -f 2 GTEx.blood.variants.inOrder.1col.annot.tsv | sed '1d') > GTEx_ge_blood.all.unsorted.tsv


sort -k2,2n -k3,3n Cartagene_ge_blood.all.unsorted.tsv > Cartagene_ge_blood.all.tsv
sort -k2,2n -k3,3n GTEx_ge_blood.all.unsorted.tsv > GTEx_ge_blood.all.tsv
bgzip Cartagene_ge_blood.all.tsv
bgzip GTEx_ge_blood.all.tsv


#remove those tmp files
rm GTEx.blood.variants.inOrder.tsv CaG.variants.inOrder.tsv


#create backups...
cp Cartagene_ge_blood.all.tsv.gz Cartagene_ge_blood.all.tsv.bkp.gz
cp GTEx_ge_blood.all.tsv.gz GTEx_ge_blood.all.tsv.bkp.gz


################################################
#GENERATE TPM MEDIAN VALUES PER DATASET :
################################################

gunzip -c Cartagene_ge_blood.all.tsv.gz | cut -f 1 > Cartagene_ge_blood.ensIDs.ordered.txt
gunzip -c GTEx_ge_blood.all.tsv.gz | cut -f 1 > GTEx_ge_blood.ensIDs.ordered.txt


#ADD THE TPM TO THE FILES AND KEEP THE ORIGINAL ORDER :

#Get the right TPM for GTEx
sed '1,2d' /lustre06/project/6065672/shared/GTex/v8/Gene_Expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct > GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.tsv
ln -s /lustre06/project/6065672/shared/Cartagene/flagship_paper/eQTLs/work_on_911_SampleFilters_after/tensorQTL/output/gwas_covariates/GTEx_v8_Expression/Matrix.504samples.RNAseq.kallisto.genes.unstranded.medianValues.withNames.txt


#DO A LEFT JOIN WITH R :
echo "library(dplyr)
#CAG
cag=read.table(\"Cartagene_ge_blood.ensIDs.ordered.txt\",h=F)
cag_tpm=read.table(\"Matrix.504samples.RNAseq.kallisto.genes.unstranded.medianValues.withNames.txt\",h=T)

names(cag)=c(\"GeneID\")
names(cag_tpm)=c(\"GeneID\",\"MedianTPM\")


#GTEX
gtex=read.table(\"GTEx_ge_blood.ensIDs.ordered.txt\",h=F)
gtex_tpm=read.table(\"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.tsv\",h=T,sep=\"\t\")

names(gtex)=c(\"GeneID\")
#mod tpm gtex
gtex_tpm=gtex_tpm[,c(\"Name\",\"Whole.Blood\")]
names(gtex_tpm)=c(\"GeneID\",\"MedianTPM\")
gtex_tpm\$GeneID=gsub(\"\\\..*\",\"\",gtex_tpm\$GeneID,perl=T)

#SOME ARE DUPLICATED, KEEP THE ONE WITH MORE TPM
gtex_tpm=gtex_tpm[order(gtex_tpm\$MedianTPM,decreasing=T),]
gtex_tpm=distinct(gtex_tpm,GeneID,.keep_all=1)

merged_cag=left_join(cag,cag_tpm,by=\"GeneID\")
write.table(merged_cag,\"Cartagene_ge_blood.ensIDs.ordered.withTPM.txt\",quote=F,sep=\"\t\",row.names=F)


merged_gtex=left_join(gtex,gtex_tpm,by=\"GeneID\")
write.table(merged_gtex,\"GTEx_ge_blood.ensIDs.ordered.withTPM.txt\",quote=F,sep=\"\t\",row.names=F)

" > addTPM.Rscript

R CMD BATCH addTPM.Rscript


#PASTE THE TPM TO CARTAGENE, THEN TO GTEx
gunzip Cartagene_ge_blood.all.tsv.gz
mv Cartagene_ge_blood.all.tsv Cartagene_ge_blood.all.tsv.bkp

paste <(cut -f 1-17 Cartagene_ge_blood.all.tsv.bkp) <(cut -f 2 Cartagene_ge_blood.ensIDs.ordered.withTPM.txt | sed '1d') <(cut -f 19 Cartagene_ge_blood.all.tsv.bkp) > Cartagene_ge_blood.all.tsv
bgzip Cartagene_ge_blood.all.tsv
tabix -s 2 -b 3 -e 3 Cartagene_ge_blood.all.tsv.gz


gunzip GTEx_ge_blood.all.tsv.gz
mv GTEx_ge_blood.all.tsv GTEx_ge_blood.all.tsv.bkp
paste <(cut -f 1-17 GTEx_ge_blood.all.tsv.bkp) <(cut -f 2 GTEx_ge_blood.ensIDs.ordered.withTPM.txt | sed '1d') <(cut -f 19 GTEx_ge_blood.all.tsv.bkp) > GTEx_ge_blood.all.tsv
bgzip GTEx_ge_blood.all.tsv
tabix -s 2 -b 3 -e 3 GTEx_ge_blood.all.tsv.gz




#next step... go back to root of fivex
cd ../../


#Generate chromosome sizes, since the browser wants to have chunks of 1Mb, and we don't want to go over the chromosome sizes, can be placed elswhere in the script...
#GENERATE INTERVALS TO EXTRACT :

echo "1 1 248956422
10 1 133797422
11 1 135086622
12 1 133275309
13 1 114364328
14 1 107043718
15 1 101991189
16 1 90338345
17 1 83257441
18 1 80373285
19 1 58617616
2 1 242193529
20 1 64444167
21 1 46709983
22 1 50818468
3 1 198295559
4 1 190214555
5 1 181538259
6 1 170805979
7 1 159345973
8 1 145138636
9 1 138394717" > chr.sizes.txt


for i in $(cat chr.sizes.txt)
do IFS=$'\n'
chr=$(echo $i | cut -d " " -f 1)
start=$(echo $i | cut -d " " -f 2)
end=$(echo $i | cut -d " " -f 3)
last=0

for f in $(seq $start 1000000 $end)
do
if [[ $last -eq 0 ]] 
then 
>&2 echo "#beginning"
elif [[ $f-$last -lt 1000000 ]] 
then 
echo "$chr $last $end" >> allIntervals.genomeWide.GRCh38.txt
else
echo "$chr $last $(($f-1))" >> allIntervals.genomeWide.GRCh38.txt
fi
last=$f
done
echo "$chr $last $end" >> allIntervals.genomeWide.GRCh38.txt
done




#NEXT STEP :
#Create the credible sets files, but without the fine mapping information for now... #TODO

cd data/cartagene_ge/

#COPY TO ebi_original dataset
cp GTEx_ge_blood.all.tsv.gz* ../ebi_original/ge/GTEx/
cp Cartagene_ge_blood.all.tsv.gz* ../ebi_original/ge/Cartagene/




#################################
######### CREDIBLE SETS #########
#################################

#GENERATE credible sets based on eQTL results as well

#phenotype_id	variant_id	chr	pos 	ref 	alt 	cs_id	cs_index	finemapped_region	pip 	z 	cs_min_r2	cs_avg_r2	cs_size	posterior_mean	posterior_sd	cs_log10bf

mkdir <path_to_eqtl_browser>/fivex/data/credible_sets/ge/Cartagene
mkdir <path_to_eqtl_browser>/fivex/data/credible_sets/ge/GTEx


#GENERATE A FIRST VERSION OF IT, WILL NEED TO HAVE 17 COLUMNS IN ORDER TO WORK...
#Cartagene
#extract the right column
gunzip -c Cartagene_ge_blood.all.tsv.gz | awk '{print $1"\t"$6"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t0\tchr"$2":"$3"-"$3"\t0\t0\t0\t0\t0\t0\t0\t0\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$18}' > <path_to_eqtl_browser>/fivex/data/credible_sets/ge/Cartagene/Cartagene.blood_ge.purity_filtered.sorted.txt
bgzip <path_to_eqtl_browser>/fivex/data/credible_sets/ge/Cartagene/Cartagene.blood_ge.purity_filtered.sorted.txt
tabix -s 3 -b 4 -e 4 <path_to_eqtl_browser>/fivex/data/credible_sets/ge/Cartagene/Cartagene.blood_ge.purity_filtered.sorted.txt.gz

#GTEx
#extract the right column
gunzip -c GTEx_ge_blood.all.tsv.gz | awk '{print $1"\t"$6"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t0\tchr"$2":"$3"-"$3"\t0\t0\t0\t0\t0\t0\t0\t0\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$18}'  > <path_to_eqtl_browser>/fivex/data/credible_sets/ge/GTEx/GTEx.blood_ge.purity_filtered.sorted.txt
bgzip <path_to_eqtl_browser>/fivex/data/credible_sets/ge/GTEx/GTEx.blood_ge.purity_filtered.sorted.txt
tabix -s 3 -b 4 -e 4 <path_to_eqtl_browser>/fivex/data/credible_sets/ge/GTEx/GTEx.blood_ge.purity_filtered.sorted.txt.gz



#REMOVE SOME COLUMNS, NEED 17 TO WORK WITH THE BROWSER (can be optimized)
cd ../credible_sets/ge/
mkdir test
cp GTEx/GTEx.blood_ge.purity_filtered.sorted.txt.gz* Cartagene/Cartagene.blood_ge.purity_filtered.sorted.txt.gz* test/



cd GTEx/ 
zcat GTEx.blood_ge.purity_filtered.sorted.txt.gz | cut -f 1-17 > GTEx.blood_ge.purity_filtered.sorted.txt
rm GTEx.blood_ge.purity_filtered.sorted.txt.gz GTEx.blood_ge.purity_filtered.sorted.txt.gz.tbi
bgzip GTEx.blood_ge.purity_filtered.sorted.txt
tabix -s 3 -b 4 -e 4 GTEx.blood_ge.purity_filtered.sorted.txt.gz 

cd ../Cartagene
zcat Cartagene.blood_ge.purity_filtered.sorted.txt.gz | cut -f 1-17 > Cartagene.blood_ge.purity_filtered.sorted.txt
rm Cartagene.blood_ge.purity_filtered.sorted.txt.gz Cartagene.blood_ge.purity_filtered.sorted.txt.gz.tbi
bgzip Cartagene.blood_ge.purity_filtered.sorted.txt
tabix -s 3 -b 4 -e 4 Cartagene.blood_ge.purity_filtered.sorted.txt.gz 

cd ../../../




#CREATE THE INDEX.TSV FILE AND COPY THE FILES AT THE RIGHT PLACE
echo "Cartagene blood credible_sets/ge/test/Cartagene.blood_ge.purity_filtered.sorted.txt.gz
GTEx blood credible_sets/ge/test/GTEx.blood_ge.purity_filtered.sorted.txt.gz" > all_EBI_credible_sets_data_index.tsv


#CREATE THE CREDIBLE SET COMBINING ALL THE DATASET, PER CHROMOSOME
for i in $(cat ../chr.sizes.txt)
do IFS=$'\n'
chr=$(echo $i | cut -d ' ' -f 1)
start=$(echo $i | cut -d ' ' -f 2)
end=$(echo $i | cut -d ' ' -f 3)
python3 ../util/merge.files.with.sorted.positions.py all_EBI_credible_sets_data_index.tsv $chr $start $end credible_sets/ge/chr$chr.ge.credible_set.tsv.gz 3
tabix -s 5 -b 6 -e 6 credible_sets/ge/chr$chr.ge.credible_set.tsv.gz
done





#dataset, tissue, trait, vchr, vpos, vref, valt, vid, ma_samples, maf, pvalue, beta, se, vtype, ac, an, r2, mol_trait_obj_id, gid, median_tpm, rsid
#BLUEPRINT       monocyte        ENSG00000116266 chr1_108776183_C_T      1       108776183       C       T       ENSG00000116266_L2      L2      chr1:107746674-109746674        0.00657144529019305     -8.10775296404622       0.665032210449112       0.986798706831566       127     -0.00348390557142348    0.0434044771475337      7.66153803775023        75      0.21466 6.0399e-13      -0.101223       0.0130229       SNP     300     382     NA      ENSG00000116266 ENSG00000116266 189.226 rs1759479       STXBP3

#at this step : 
#Cartagene       blood   ENSG00000227232 chr1_16949_A_C  1       16949   A       C       chr1_16949_A_C  0       chr1:16949-16949        0       0       0       0       0       0       0       0
#would need more : 
#0	0	0	0	0	TYPE	0	0	NA	$3	$3	0	RS	GENEID

cd credible_sets/ge/
mkdir temp
mv chr*tsv.gz* temp
cd temp

#with R :
#get gene ids and rs ids to annotate the files : 

for f in $(seq 1 22)
do
zcat chr$f.ge.credible_set.tsv.gz | cut -f 3,4 > chr$f.ensIDs_Pos.txt 

echo "library(dplyr)
test=read.table(\"chr$f.ensIDs_Pos.txt\",h=F)
names(test)=c(\"ENSID\",\"chr_pos_ref_alt\")
gene_names=read.table(\"<path_to_eqtl_files>/geneConversion.txt\",h=F)
annotations=read.table(\"<path_to_eqtl_browser>/fivex/data/cartagene_ge/CaG.allVariants.uniq.dbSNP.annot.2cols.tsv\",h=F)
names(annotations)=c(\"chr_pos_ref_alt\",\"id\")
names(gene_names)=c(\"ENSID\",\"HUGO\")

annotations\$chr_pos_ref_alt=gsub(\"^\",\"chr\",annotations\$chr_pos_ref_alt)

merged_test_names=left_join(test,gene_names,by=\"ENSID\")
merged_test_names_rs=left_join(merged_test_names,distinct(annotations,chr_pos_ref_alt,.keep_all=1),by=\"chr_pos_ref_alt\")
write.table(merged_test_names_rs,\"chr$f.ensIDs_Pos.withSymbol_rs.txt\")" > chr$f.Rscript

R CMD BATCH chr$f.Rscript
done





#SHOULD LOOK LIKE (33 columns file): 

#1) study: str
#2) tissue: str
## The rest of these fields are present in all credible interval files
#3) gene_id: str  # this column is labeled "phenotype_id" in the original file
#4) var_id: str  # in chrom_pos_ref_alt format -- not used
#5) chromosome: str
#6) position: int
#7) ref_allele: str
#8) alt_allele: str
#9) cs_id: str
#10) cs_index: str
#11) finemapped_region: str
#12) pip: float
#13) z: float
#14) cs_min_r2: float
#15) cs_avg_r2: float
#16) cs_size: int
#17) posterior_mean: float
#18) posterior_sd: float
#19) cs_log10bf: float
## Extra fields after data joining
#20) ma_samples: ty.Optional[int] = None
#21) maf: ty.Optional[float] = None
#22) log_pvalue: ty.Optional[float] = None
#23) beta: ty.Optional[float] = None
#24) stderr_beta: ty.Optional[float] = None
#25) type: ty.Optional[str] = None
#26) ac: ty.Optional[int] = None
#27) an: ty.Optional[int] = None
## The r2 field may contain 'NA's -- see note in VariantContainers
#28) r2: ty.Optional[float] = None
#29) mol_trait_obj_id: ty.Optional[str] = None
#30) gid: ty.Optional[str] = None
#31) median_tpm: ty.Optional[float] = None
#32) rsid: ty.Optional[str] = None
#33) symbol: ty.Optional[str] = None


# ADD THE VARIANT TYPE NOW...
for f in $(seq 1 22)
do
paste <(zcat chr$f.ge.credible_set.tsv.gz | cut -f 1-19 | awk '{print $0}') <(zcat chr$f.ge.credible_set.tsv.gz | awk '{print $20"\t"$21"\t"$22"\t"$23"\t"$24}') <(zcat chr$f.ge.credible_set.tsv.gz | awk '{ type="SNP" ; if(length($7) > 1 || length($8) > 1 ) type="INDEL" ; print type}' ) <(zcat chr$f.ge.credible_set.tsv.gz | awk '{print "0\t0\tNA\t"$3"\t"$3"\t"$25}') <(sed -e '1d' -e 's/"//g' -e 's/ /\t/g' chr$f.ensIDs_Pos.withSymbol_rs.txt | awk '{print $5"\t"$4}') > ../chr$f.ge.credible_set.tsv
done

cd ../
for f in chr*tsv ; do bgzip $f  ; tabix -s 5 -b 6 -e 6 $f.gz ; done

#What we have (33 columns)... 
#Cartagene       blood   ENSG00000223972 chr1_13649_G_C  1       13649   G       C       chr1_13649_G_C  0       chr1:13649-13649        0       0       0       0       0       0       0       0       0       0           0       0       0       SNP     0       0       NA      ENSG00000223972 ENSG00000223972 0	rs11452411	SCYL3




#DO NEXT STEP THE ebi_ge folder... per chunk of 1MB :
cd ../../../



for i in $(cat allIntervals.genomeWide.GRCh38.txt)
do IFS=$'\n'
chr=$(echo $i | cut -d ' ' -f 1)
start=$(echo $i | cut -d ' ' -f 2)
end=$(echo $i | cut -d ' ' -f 3)
python util/merge.files.with.sorted.positions.py index_projects.txt  $chr $start $end data/ebi_ge/$chr/all.EBI.ge.data.chr$chr.$start-$end.tsv.gz 2
done


#module load tabix
#do the index files after... otherwise there's an issue because the machine is too fast to generate them.
cd <path_to_eqtl_browser>/fivex/data/ebi_ge

for i in $(seq 1 22)
do cd $i
for f in *.tsv.gz
do
tabix -f -s 4 -b 5 -e 5 $f
done
cd <path_to_eqtl_browser>/fivex/data/ebi_ge
done





#CREATE FINAL DATABASE :
#SQLITE3 database :

#example : table sig
#header seems to be : 
#pvalue|dataset|tissue|GENE|chr|pos|ref|alt|NA|CS|CSsize
#0.00657144529019305|BLUEPRINT_SE|monocyte|ENSG00000116266|1|108776183|C|T|L2|127
#0.0150511076419484|Quach_2016|monocyte_R848|ENSG00000116266|1|108777306|A|G|L1|18
#0.00657144529019305|BLUEPRINT_SE|monocyte|ENSG00000116266|1|108777768|T|C|L2|127
#0.00657144529019305|BLUEPRINT_SE|monocyte|ENSG00000116266|1|108778306|T|G|L2|127
#0.00657144529019305|BLUEPRINT_SE|monocyte|ENSG00000116266|1|108778623|G|A|L2|127
#0.00872929319149929|BLUEPRINT_SE|monocyte|ENSG00000116266|1|108778690|T|TTATA|L2|127
#0.00657144529019305|BLUEPRINT_SE|monocyte|ENSG00000116266|1|108778722|A|T|L2|127
#0.00657144529019305|BLUEPRINT_SE|monocyte|ENSG00000116266|1|108778735|T|C|L2|127


#ENSG00000227232 1       10433   ACCCTAAC        ACCCCTAAC       chr1_10433_ACCCTAAC_ACCCCTAAC   56      0.0575396865606308      0.00253663452912465     0.403804987668991    0.132967531681061       INDEL   58      1007    1       ENSG00000227232 ENSG00000227232 0       chr1_10433_ACCCTAAC_ACCCCTAAC

cd <path_to_eqtl_browser>/fivex/data/cartagene_ge

#temporary files
mkdir for_sql
gunzip -c Cartagene_ge_blood.all.tsv.gz | awk '{print "-\tCartagene\tblood\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t0\t0\t"$9}' > for_sql/Cartagene_ge_blood.all.forSQL.tsv

gunzip -c GTEx_ge_blood.all.tsv.gz | awk '{print "-\tGTEx\tblood\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t0\t0\t"$9}' > for_sql/GTEx_ge_blood.all.forSQL.tsv
cd for_sql

#combine all the datasets
cat Cartagene_ge_blood.all.forSQL.tsv GTEx_ge_blood.all.forSQL.tsv | sort -k5,5n -k6,6n > forSQL.all.tsv


#DO THE DATABASE :
sqlite3 pip.best.variant.summary.sorted.indexed.sqlite3.db

#> PRAGMA index_list(sig);
#0|idx_gene|0|c|0
#1|idx_study_tissue|0|c|0
#2|idx_chrom_pos|0|c|0

create table sig(pip real, study text, tissue text, gene_id text, chrom text, pos int, ref text, alt text, cs_index text, cs_size int, pvalue real);
.mode tabs
.import forSQL.all.tsv sig

CREATE INDEX idx_gene ON sig (gene_id);
CREATE INDEX idx_study_tissue ON sig (study, tissue);
CREATE INDEX idx_chrom_pos ON sig (chrom, pos);
#exit

#copy at the right place!
cp pip.best.variant.summary.sorted.indexed.sqlite3.db ../../credible_sets/ge/pip.best.variant.summary.sorted.indexed.sqlite3.db




