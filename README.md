#RNAseq results treatment and eQTL calling script

#1) Process RNA-seq
bash 01_RNASeq_processing_unstranded_kallisto.sh

#2) Form the Count table using tximport
Rscript 02_kalliso_tximport.R

#3) Filter RNAseq, batch normalize, remove mismatches for sex and INT normalize

Rscript 03_RNAseq_filtering_and_normalization.R

#4) Run PEER, to correct for unwanted variations
#Run the R session from a singularity container
module load StdEnv/2020
module load apptainer

singularity run --bind <path_to_data>:/data <path_to_container>/peer_data.sif Rscript 04_PEER_correction_step1.R

#5) Look for effect of PEER correction
Rscript 05_Look_for_outliers_post_PEER.R

#6) Apply a second round of PEER correction, without the outliers detected
singularity run --bind <path_to_data>:/data <path_to_container>/peer_data.sif Rscript 06_PEER_correction_step2.R

#7) Verify PEER results
Rscript 07_Verify_PEER.R

#Remove outlier samples and samples with sex mismatches from the original matrix and from the WGS files, ids will be in samplesToRemove.txt file
grep -Fwvf samplesToRemove.txt RNAseq.files.runs.911.txt > RNAseq.files.runs.911.noBadSamples.txt

#8) Extract common samples from RNAseq and WGS and prepare the files for eQTL analysis
bash 08_Extract_commonSamples_RNAseq_WGS.sh

#9) Launch eQTL analysis with tensorQTL
bash 09_Launch_tensorQTL.sh

#From this point, multiple processing of files are needed.
#10) Liftover positions from GTEx v9 results to b38 coordinates
bash 10_liftover_gtex_b37_to_b38.sh

#11) Parse results from tensorQTL and compare them to GTEx results, eQTL-gen results, and form the input for the database
Rscript 11_compare_hits_to_gtex.allgenes.504.03112024_fdr_per_gene.withHLA.R

#12) Parse for the eQTL fiveX database

bash 12_Parse_results_eQTL_browser.sh


