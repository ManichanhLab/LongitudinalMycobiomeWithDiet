#! usr/env/bash     

sampdir="/mnt/synology/ALEIX/c_albicans_inst_pasteur/Mapping"
outdir="/mnt/synology/ALEIX/c_albicans_inst_pasteur/Mapping"
scriptdir="/mnt/storage5TB/zixuan"
mgdir="/mnt/synology/DATABASES/FungiDB/4Kmgs_clstr"
mgeach="/mnt/synology/ZIXUAN/GENOMIC/BUSCO/Candida/busco"
pdir="/mnt/synology/DATABASES/QIIME2/ITS/Unite_Refseq"
pindex="/mnt/synology/ALEIX"

samplefile=$1     # This is the file which contains the list of samples to process.         

while read p 
do 
#mkdir $outdir/${p}/${p}
#mkdir $outdir

#echo "Start mapping ${p} on 4Kmgs"

bowtie2-build -f $outdir/${p}/${p}"_ITS.fasta" $outdir/${p}/${p}"index"
bowtie2-build -f $mgeach/${p}/"busco_sequences"/${p}"_all_single.fna" $mgeach/${p}/"busco_sequences"/${p}"index"

bowtie2 -p 40 -x $mgeach/${p}/"busco_sequences"/${p}"index" \
-1 $sampdir/${p}/"reads_paired_1.fastq" \
-2 $sampdir/${p}/"reads_paired_2.fastq" > $outdir/${p}/${p}_MG.sam

bowtie2 -p 40 -x $outdir/${p}/${p}"index" \
-1 $sampdir/${p}/"reads_paired_1.fastq" \
-2 $sampdir/${p}/"reads_paired_2.fastq" > $outdir/${p}/${p}_ITS.sam

#bwa mem $mgdir/"clstr_4Kmgs_clstr.fasta" \
#$sampdir/${p}"_paired_1.fastq" $sampdir/${p}"_paired_2.fastq" > $outdir/${p}/${p}.sam  
#bwa mem $sampdir/"allFungi_seqs_qiime.fasta" \
#$sampdir/${p}"_paired_1.fastq" $sampdir/${p}"_paired_2.fastq" > $outdir/${p}/${p}.sam



# bowtie2 -p 30 -x $mgdir/clstr_4Kmgs \
# -q $sampdir/${p}_kneaddata.fastq.gz \
# -S $outdir/${p}/${p}.sam 2> $outdir/${p}/${p}.log

samtools sort -n -O sam $outdir/${p}/${p}_MG.sam | samtools fixmate -m -O bam - $outdir/${p}/${p}_MG.bam
rm $outdir/${p}/${p}_MG.sam
samtools sort -O bam -o $outdir/${p}/${p}_MG.sorted.bam $outdir/${p}/${p}_MG.bam
rm $outdir/${p}/${p}_MG.bam
samtools view -h -b -q 30 $outdir/${p}/${p}_MG.sorted.bam > $outdir/${p}/${p}_MG.sorted_q30.bam
samtools sort -n -O sam $outdir/${p}/${p}_ITS.sam | samtools fixmate -m -O bam - $outdir/${p}/${p}_ITS.bam
rm $outdir/${p}/${p}_ITS.sam
samtools sort -O bam -o $outdir/${p}/${p}_ITS.sorted.bam $outdir/${p}/${p}_ITS.bam
rm $outdir/${p}/${p}_ITS.bam
samtools view -h -b -q 30 $outdir/${p}/${p}_ITS.sorted.bam > $outdir/${p}/${p}_ITS.sorted_q30.bam

#samtools view -Sq 30 $outdir/${p}/${p}.sam > $outdir/${p}/${p}.30.sam
#samtools view -bSq 30 $outdir/${p}/${p}.sam > $outdir/${p}/${p}.30.bam
#samtools sort -o $outdir/${p}/${p}.30.sorted.bam $outdir/${p}/${p}.30.bam
#python3 $scriptdir/coverageFilter.py -i $outdir/${p}/${p}.30.sam -o $outdir/${p}/${p}.30.c80.list
#samtools view $outdir/${p}/${p}.30.sorted.bam | grep -f $outdir/${p}/${p}.30.c80.list > $outdir/${p}/${p}.30.c80.sam 
#samtools view -bt $mgdir/4Kmgs_clstr.fasta.fai -o $outdir/${p}/${p}.30.c80.bam $outdir/${p}/${p}.30.c80.sam 
#samtools index $outdir/${p}/${p}.30.c80.bam
#samtools idxstats $outdir/${p}/${p}.30.c80.bam > $outdir/${p}/${p}.idxstats.txt
samtools depth $outdir/${p}/${p}_MG.sorted_q30.bam | gzip > $outdir/${p}/${p}_MG.depth.30.txt.gz

#samtools view -Sq 30 $outdir/${p}/${p}.sam > $outdir/${p}/${p}.30.sam
#samtools view -bSq 30 $outdir/${p}/${p}.sam > $outdir/${p}/${p}.30.bam
#samtools sort -o $outdir/${p}/${p}.30.sorted.bam $outdir/${p}/${p}.30.bam
#python3 $scriptdir/coverageFilter.py -i $outdir/${p}/${p}.30.sam -o $outdir/${p}/${p}.30.c80.list
#samtools view $outdir/${p}/${p}.30.sorted.bam | grep -f $outdir/${p}/${p}.30.c80.list > $outdir/${p}/${p}.30.c80.sam 
#samtools view -bt $outdir/${p}/allFungi_seqs_qiime.fasta.fai -o $outdir/${p}/${p}.30.c80.bam $outdir/${p}/${p}.30.c80.sam 
#samtools index $outdir/${p}/${p}.30.c80.bam
#samtools idxstats $outdir/${p}/${p}.30.c80.bam > $outdir/${p}/${p}.idxstats.txt
samtools depth $outdir/${p}/${p}_ITS.sorted_q30.bam | gzip > $outdir/${p}/${p}_ITS.depth.30.txt.gz 

#    if [[ -f "$outdir/${p}/${p}.idxstats.txt" ]]
 #   then   # When the pipeline has been executed properly and the script recognizes the output-
  #       printf "removing secondary files files for ${p}\n\n"
 #        #rm $outdir/${p}/${p}.30.sorted.bam
  #       rm $outdir/${p}/${p}.sam
 #   else echo "No idxstats generated"
 #   fi


#gawk -i inplace '$3 > 0' $outdir/${p}/${p}".idxstats.txt"


#sed 's/|/\t/g' $outdir/${p}/${p}".idxstats.txt" | awk '{print $1,$4}' > $outdir/${p}/${p}"_assembly_hp.txt"
#python3 $scriptdir/make_buglist.py -i $outdir/${p}/${p}"_assembly_hp.txt" -o $outdir/${p}/${p}"_buglist.txt"
#python3 $scriptdir/buglist.py -i $outdir/${p}/${p}"_assembly_hp.txt" -o $outdir/${p}/${p}"_buglist_stratified.txt"
#gawk -i inplace -F'\t' 'BEGIN{OFS="\t"}{print $2,$3}' $outdir/${p}/${p}"_buglist_stratified.txt"

done < $samplefile
