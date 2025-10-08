#liftover gtex
gunzip -c Whole_Blood.allpairs.txt.gz | cut -f 2 > Whole_Blood.allpairs.pos.bed

sed -i 's/^/chr/' Whole_Blood.allpairs.pos.bed
paste <(sed -e 's/_/\t/g' -e 's/\tb37//' Whole_Blood.allpairs.pos.bed | awk '{print $1"\t"$2-1"\t"$2}') <(cat Whole_Blood.allpairs.pos.bed) > Whole_Blood.allpairs.pos.hg19.bed
sed -i '1d' Whole_Blood.allpairs.pos.hg19.bed

liftOver Whole_Blood.allpairs.pos.hg19.bed ~/projects/def-hussinju/shared/References/liftOver/hg19ToHg38.over.chain Whole_Blood.allpairs.pos.GRCh38.bed Whole_Blood.allpairs.pos.hg19.unmapped


paste <(cut -f 1,3 Whole_Blood.allpairs.pos.GRCh38.bed) <(cut -f 4 Whole_Blood.allpairs.pos.GRCh38.bed | sed -e 's/_/\t/' -e 's/_/\t/' -e 's/_b37//'| cut -f 3) | sed 's/\t/_/g' | uniq > Whole_Blood.allpairs.GRCh38.newIDs.txt

paste <(cut -f 4 Whole_Blood.allpairs.pos.GRCh38.bed) Whole_Blood.allpairs.GRCh38.newIDs.txt | sed 's/^chr//' > Whole_Blood.allpairs.GRCh38.old_vs_newIDs.txt

#merge new positions to the original file
R

library(data.table)

library(dplyr)
b37=as.data.frame(fread("Whole_Blood.allpairs.txt.gz",h=T))
b38_ids=as.data.frame(fread("Whole_Blood.allpairs.GRCh38.old_vs_newIDs.sorted.txt",h=F))
names(b38_ids)=c("variant_id","variant_id_b38")

b37_necessary=b37[,c("gene_id","variant_id","maf","pval_nominal","slope")]
rm(b37)

merged_gtex=left_join(b37_necessary,b38_ids,by="variant_id")

write.table(merged_gtex,file="Whole_Blood.allpairs.withB38_ids.txt")


