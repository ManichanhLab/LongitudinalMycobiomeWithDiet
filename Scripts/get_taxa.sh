
outdir="/mnt/synology/ALEIX/"

genusname=$1


while IFS= read -r -u3 line
do
    
   # esearch -db nuccore -query "$line" | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS="\t";}{print $3;}' 
   # esearch -db taxonomy -query "$line" | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS="\t";}{print $3;}' 
    esearch -db SRA -query $line | elink -target taxonomy | efetch -format native -mode xml | grep ScientificName | awk -F ">|<" 'BEGIN{ORS="\t";}{print $3;}'
    echo $line 
done 3< $genusname >> $outdir/calb_taxa_full.txt

