#! usr/env/bash  
#run in /mnt/synology/busco_mapping/hmm/ folder
#INPUT: FULL NAME OF GENOMES TO BE ANALYZED (including .fasta/.fa/...)
hmm_dir="/mnt/synology/ALEIX/MASTER/ITS"
genomes_dir="/mnt/synology/ALEIX/FUNGAL_GENOMES/NCBI"

genomes=$1

while read p
do
	mkdir "$p"
	nhmmer --noali -A  "$p"/"$p"_final.aln --tblout "$p"/"$p"_final.tb $hmm_dir/final.hmm $genomes_dir/"$p" > "$p"/"$p"_output_final
	awk '($1!~"#" && $13 < 0.001){print $1,$3,$7,$8,$12,$13}' "$p"/*_final.tb | sort -k1  > "$p"/"$p"_filtered_output_final.tbl
done < $genomes
python3 /mnt/synology/ALEIX/scripts/cnv.py $genomes