#! usr/env/bash     

sampdir="/mnt/synology/ALEIX/mock_genome"
outdir="/mnt/synology/ALEIX/mock_genome"
scriptdir="/mnt/storage5TB/zixuan"
mgdir="/mnt/synology/DATABASES/FungiDB/4Kmgs_clstr"
mgeach="/mnt/synology/ZIXUAN/GENOMIC/BUSCO/NCBI"
pdir="/mnt/synology/DATABASES/QIIME2/p/Unite_Refseq"
pindex="/mnt/synology/ALEIX"

p=$1     # This is the file which contains the list of samples to process.         

#while read p 
#do 
#mkdir $outdir/${p}
#mkdir $outdir

echo "Start mapping ${p} on 4Kmgs"
bowtie2-build -f $outdir/${p}"_ITS.fasta" $outdir/${p}"index"
bowtie2-build -f $mgeach/${p}/"busco_sequences"/${p}"_all_single.fna" $mgeach/${p}/"busco_sequences"/${p}"index"

bowtie2 -p 40 -x $mgeach/${p}/"busco_sequences"/${p}"index" \
-1 $sampdir/"yjm693_mock_paired_1.fastq" \
-2 $sampdir/"yjm693_mock_paired_2.fastq" > $outdir/${p}_MG.sam

bowtie2 -p 40 -x $outdir/${p}"index" \
-1 $sampdir/"yjm693_mock_paired_1.fastq" \
-2 $sampdir/"yjm693_mock_paired_2.fastq" > $outdir/${p}_ITS.sam

#bwa mem $mgdir/"clstr_4Kmgs_clstr.fasta" \
#$sampdir/${p}"_paired_1.fastq" $sampdir/${p}"_paired_2.fastq" > $outdir/${p}.sam  
#bwa mem $sampdir/"allFungi_seqs_qiime.fasta" \
#$sampdir/${p}"_paired_1.fastq" $sampdir/${p}"_paired_2.fastq" > $outdir/${p}.sam



# bowtie2 -p 30 -x $mgdir/clstr_4Kmgs \
# -q $sampdir/${p}_kneaddata.fastq.gz \
# -S $outdir/${p}.sam 2> $outdir/${p}.log

samtools sort -n -O sam $outdir/${p}_ITS.sam | samtools fixmate -m -O bam - $outdir/${p}_ITS.bam
rm $outdir/${p}_ITS.sam
samtools sort -O bam -o $outdir/${p}_ITS.sorted.bam $outdir/${p}_ITS.bam
rm $outdir/${p}_ITS.bam
samtools view -h -b -q 30 $outdir/${p}_ITS.sorted.bam > $outdir/${p}_ITS.sorted_q30.bam
samtools sort -n -O sam $outdir/${p}_MG.sam | samtools fixmate -m -O bam - $outdir/${p}_MG.bam
rm $outdir/${p}_MG.sam
samtools sort -O bam -o $outdir/${p}_MG.sorted.bam $outdir/${p}_MG.bam
rm $outdir/${p}_MG.bam
samtools view -h -b -q 30 $outdir/${p}_MG.sorted.bam > $outdir/${p}_MG.sorted_q30.bam

#samtools view -Sq 30 $outdir/${p}_MG.sam > $outdir/${p}_MG.30.sam
#samtools view -bSq 30 $outdir/${p}_MG.sam > $outdir/${p}_MG.30.bam
#samtools sort -o $outdir/${p}_MG.30.sorted.bam $outdir/${p}_MG.30.bam
#python3 $scriptdir/coverageFilter.py -i $outdir/${p}_MG.30.sam -o $outdir/${p}_MG.30.c80.list
#samtools view $outdir/${p}_MG.30.sorted.bam | grep -f $outdir/${p}_MG.30.c80.list > $outdir/${p}_MG.30.c80.sam 
#samtools view -bt $mgdir/4Kmgs_clstr.fasta.fai -o $outdir/${p}_MG.30.c80.bam $outdir/${p}_MG.30.c80.sam 
#samtools index $outdir/${p}.30.c80.bam
#samtools idxstats $outdir/${p}.30.c80.bam > $outdir/${p}.idxstats.txt
samtools depth $outdir/${p}_MG.sorted_q30.bam | gzip > $outdir/${p}_MG.depth.30.txt.gz

#samtools view -Sq 30 $outdir/${p}_ITS.sam > $outdir/${p}_ITS.30.sam
#samtools view -bSq 30 $outdir/${p}_ITS.sam > $outdir/${p}_ITS.30.bam
#samtools sort -o $outdir/${p}_ITS.30.sorted.bam $outdir/${p}_ITS.30.bam
#python3 $scriptdir/coverageFilter.py -i $outdir/${p}_ITS.30.sam -o $outdir/${p}_ITS.30.c80.list
#samtools view $outdir/${p}_ITS.30.sorted.bam | grep -f $outdir/${p}_ITS.30.c80.list > $outdir/${p}_ITS.30.c80.sam 
#samtools view -bt $outdir/${p}/allFungi_seqs_qiime.fasta.fai -o $outdir/${p}_ITS.30.c80.bam $outdir/${p}_ITS.30.c80.sam 
#samtools index $outdir/${p}.30.c80.bam
#samtools idxstats $outdir/${p}.30.c80.bam > $outdir/${p}.idxstats.txt
samtools depth $outdir/${p}_ITS.sorted_q30.bam | gzip > $outdir/${p}_ITS.depth.30.txt.gz 

#    if [[ -f "$outdir/${p}.idxstats.txt" ]]
 #   then   # When the pipeline has been executed properly and the script recognizes the output-
  #       printf "removing secondary files files for ${p}\n\n"
 #        #rm $outdir/${p}.30.sorted.bam
  #       rm $outdir/${p}.sam
 #   else echo "No idxstats generated"
 #   fi


#gawk -i inplace '$3 > 0' $outdir/${p}".idxstats.txt"


#sed 's/|/\t/g' $outdir/${p}".idxstats.txt" | awk '{print $1,$4}' > $outdir/${p}"_assembly_hp.txt"
#python3 $scriptdir/make_buglist.py -i $outdir/${p}"_assembly_hp.txt" -o $outdir/${p}"_buglist.txt"
#python3 $scriptdir/buglist.py -i $outdir/${p}"_assembly_hp.txt" -o $outdir/${p}"_buglist_stratified.txt"
#gawk -i inplace -F'\t' 'BEGIN{OFS="\t"}{print $2,$3}' $outdir/${p}"_buglist_stratified.txt"

#done < $samplefile
