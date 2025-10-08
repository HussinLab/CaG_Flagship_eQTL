#REPROCESS RNA-SEQ UNSTRANDED:

#for each run...
for run in run134 run135 run147-148-159 run151 run152 run166 run167 run185 run190 run191 run203 run204 run229 run230 run254 run255 run284 run285 run289 run290 run310 run311 run312 run313 run319 run320 run322 run323 run339 run340 run481 run510 run511 run535 run539 run540 run589
cd $run/trimmed_files

echo '#!/bin/bash' > scriptKallisto.E98.unstranded.sh
echo "module load StdEnv/2020
module load kallisto/0.46.1

file=\$1
gtf=<path_to_gtf>/Homo_sapiens.GRCh38.98.gtf
ref=<path_to_kallisto_idx_E98>/Homo_sapiens.GRCh38.cdna.all_kallisto_index
mkdir \$file.GRCh38.E98.unstranded.kallisto
kallisto quant -i \$ref -o \$file.GRCh38.E98.unstranded.kallisto -t 6 -b 100 --bias -g \$gtf \${file}_R1_val_1.fq.gz \${file}_R2_val_2.fq.gz" >> scriptKallisto.E98.unstranded.sh

#launch parallel jobs
echo '#!/bin/bash' > launchKallisto.sh
echo "parallel --jobs 8 ./scriptKallisto.E98.unstranded.sh {1} :::: files.txt" >> launchKallisto.sh

chmod ug+rwx scriptKallisto.E98.unstranded.sh 
sbatch --chdir $(pwd) --account=ctb-hussinju --time=2:00:00 --output=launchKallisto.%j.out --error=launchKallisto.%j.err --cpus-per-task=48 --mem=100g launchKallisto.sh


#3 runs needs to be concatenated at the FastQ step : 147 and 148 with the run 159. So in the end, one resulting run from this step : run147-148-159