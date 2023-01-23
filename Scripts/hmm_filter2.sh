#! usr/env/bash  
#run in /mnt/synology/busco_mapping/hmm/ folder
#INPUT: FULL NAME OF GENOMES TO BE ANALYZED (including .fasta/.fa/...)
hmm_dir="/mnt/synology/ALEIX/MASTER/ITS"
genomes_dir="/mnt/synology/ALEIX/FUNGAL_GENOMES/NCBI"

genomes=$1

while read p
do
	awk '($1!~"#" && $13 < 0.001){print $1,$3,$7,$8,$12,$13}' "$p"/*_final.tb | sort -k1  > "$p"/"$p"_filtered2_output_final.tbl
done < $genomes