#!/usr/bin/env bash

sampdir="/mnt/synology/ZIXUAN/LongitudinalMycobiome/SEQUENCE"   ## RAW READS PATH
cleandir="/mnt/synology/ZIXUAN/KNEADDATA/PSO_USA"				## KNEADDATA OUTPUT PATH
outdir="/mnt/synology/ZIXUAN/LongitudinalMycobiome/ANALYSIS"	## PROFILING OUTPUT PATH
## path for kneaddata
hg37_dir1="/mnt/synology/DATABASES/KNEADDATA/HUM_GENOME"  		## KNEADDATA DB
threads=10
processes=10
trimmomatic="/home/zixuan/Downloads/trimmomatic/Trimmomatic-0.39"	## REMEMBER TO CHANGE THE PATH TO YOUR TRIMMOMATIC
trim_sliding="SLIDINGWINDOW:4:20" 
minimum_length="MINLEN:59" 
## path for humann
metaphlan="/home/gsergom/anaconda3/bin/metaphlan"				## REMEMBER TO CHANGE THE PATH TO YOUR METAPHLAN
protein_db="/mnt/synology/DATABASES/HUMANN/uniref/uniref"
choco_db="/mnt/synology/DATABASES/HUMANN/chocophlan_v296_20191b"


samplefile=$1     # This is the file which contains the list of samples to process.         

while read p 
do 
# mkdir $outdir
 mkdir $outdir/${p}
 mkdir $outdir/${p}"/Mycobiome"
 mkdir $outdir/${p}"/Microbiome"


printf "Start running FunOMICv2\n\n"

#### The following is only needed if the raw reads for each sample are splitted, COMMENT if you do not need
  printf "Start merging reads \n"  
  zcat $sampdir/${p}/*1.fq.gz > $sampdir/${p}/${p}"_1.fq.gz"
  zcat $sampdir/${p}/*2.fq.gz > $sampdir/${p}/${p}"_2.fq.gz"


#######################################
############# KNEADDATA ###############
#######################################
printf "Starting quality filtering and decontamination of ${p}\n"

  mkdir -pv "$cleandir/${p}"  
  mkdir -pv "$cleandir/logs/"

  time kneaddata --input $sampdir/${p}/${p}"_1.fq.gz" \
    --input $sampdir/${p}/${p}"_2.fq.gz" --output-prefix "${p}_kneaddata"\
    -db $hg37_dir1 -t $threads -p $processes \
    --trimmomatic $trimmomatic \
    --trimmomatic-options $trim_sliding \
    --trimmomatic-options $minimum_length \
    --output $cleandir"/"${p} \
    --log $cleandir/"logs/"${p}".log" # command to run kneaddata


  printf "output stored in: "$cleandir"/"${p}

  # The following section is just defensive programming. It checks for the most usual problems and reports it, to make it easier to report whenever something went wrong.
  if [[ -f "$cleandir/${p}/${p}_kneaddata_paired_1.fastq" && -f "$cleandir/${p}/${p}_kneaddata_paired_2.fastq" ]]; then   # When the pipeline has been executed properly and the script recognizes the output-
    printf "removing secondary files files for ${p}\n"
    find $cleandir/${p}/ -type f -not -name "*_kneaddata_paired_1.fastq" -not -name "*_kneaddata_paired_2.fastq" | xargs -d"\n" -exec rm # removes everything excep $sample_kneaddata_paired_(1|2).fastq
    if [[ ! `find $cleandir/${p}/ -type f -not -name *paired* 2>/dev/null | wc -l` -gt 0 ]]; then # checks if the temp files are found (they should have been removed by this step)
      if [[ `find $cleandir/${p}/ -type f -name *paired* 2>/dev/null | wc -l` -gt 0 ]]; then
        printf "temp files removed and clean fastq correct, compressing clean reads...\n\n"
        gzip "$cleandir/${p}/"*_paired_1.fastq "$cleandir/${p}/"*_paired_2.fastq
      else
        printf "Post-QC clean reads found before and now not anymore. This could mean they got removed. Check your sample names\n\n"
      fi
    else
      printf "Temp files could not be removed but clean fastqs found, they will be compressed and data might be still usable, but remember to remove intermediate files by hand\n\n"
    fi
    if [[ -f "$cleandir/${p}/${p}_kneaddata_paired_1.fastq.gz" ]] && [[ -f "$cleandir/${p}/${p}_kneaddata_paired_2.fastq.gz" ]]; then # checks if one or both paired.fastqs does not exist
      printf "clean fastqs for sample ${p} compressed; moving to next sample:\n\n"
    else
      printf "one or both fastqs could not be compressed. This should not be a great problem, but beware of disk space.\n\n"
    fi
  else
    printf "post-QC output not found. Check if it is actually being produced. if it does, it may be a name incompatibility, stop the script and change it, or remove the unnecessary files by hand. Your choice.\n\n"
  fi

#######################################
############## FunOMIC ################
#######################################

printf "Start mycobiome profiling \n of ${p}"

  time FunOMIC2.sh -1 "$cleandir/${p}/${p}_kneaddata_paired_1.fastq.gz" -2 "$cleandir/${p}/${p}_kneaddata_paired_2.fastq.gz" \
  -p ${p} -o $outdir/${p}"/Mycobiome" \
  -a /mnt/synology/DATABASES/UHGG/bowtie2 \
  -b /mnt/synology/DATABASES/FungiDB/FunOMIC-T/FunOMIC-T.v2 \
  -c /mnt/synology/DATABASES/FungiDB/FunOMIC-P/FunOMIC.P.v2 -t 35 

######################################
############## HUMAnN ################
######################################

printf "Start microbiome profiling \n of ${p}"

  # Merge read pairs into single fastq
  zcat "$cleandir/${p}/"*.fastq.gz >> "$cleandir/${p}/${p}_combined_reads.fastq"
  printf "merged reads stored in: $cleandir/${p}/${p}_combined_reads.fastq \n"

  #run HUMAnN 3.0
  printf "saving output in: "$outdir"/"${p}"/"${p}"_HUMAnN_results\n"

  time humann -r -i $cleandir"/"${p}"/"${p}"_combined_reads.fastq" -o $outdir"/"${p}/Microbiome \
  --protein-database $protein_db \
  --nucleotide-database $choco_db \
  --threads 35 > "$outdir/${p}/Microbiome/${p}_HUMAnN_std_output.log"

  # Remove the merged read file when finished.
  rm -f "$cleandir/${p}/${p}_combined_reads.fastq"

  # Move all relevant final output to the same folder
  mv $outdir"/"${p}/Microbiome"/"${p}"_combined_reads_humann_temp/"${p}"_combined_reads_metaphlan_bugs_list.tsv" $outdir"/"${p}/Microbiome
  mv $outdir"/"${p}/Microbiome"/"${p}"_combined_reads_humann_temp/"${p}"_combined_reads.log" $outdir"/"${p}/Microbiome

  # check if the intermediate files are present and if so, remove them.

  if [[ -f "$outdir/$sample/${p}/Microbiome/${p}_combined_reads_metaphlan_bugs_list.tsv" ]] && [[ -f "$outdir/$sample/${p}/Microbiome/${p}_combined_reads.log" ]]; then
    printf "bugs list and log succesfully moved. Checking and removing temp folder..."
    if [[ -d "$outdir/${p}/Microbiome/${p}_combined_reads_humann_temp" ]]; then
        rm -r "$outdir/${p}/Microbiome/${p}_combined_reads_humann_temp"
        printf "temp files removed succesfully"
      else
        printf "The temp files was not detected. Please, check the script for bugs"
    fi
  else
    printf "trmp files are not found"
  fi

done < $samplefile
