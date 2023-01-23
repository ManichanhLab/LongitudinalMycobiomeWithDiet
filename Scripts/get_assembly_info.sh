#! usr/env/bash

outdir="/mnt/synology/ALEIX/"

genusname=$1


while IFS=, read -r -u3 gca
do	
    
    #esearch -db assembly -query $gca | elink -target biosample | efetch -format native -mode xml | grep -oP '(?<=Sample name">).*?(?=</Id>)'
    esearch -db assembly -query $gca | esummary | grep Sub_value | awk -F ">|<" 'BEGIN{ORS="\t";}{print $3;}'
    echo $gca
done 3< $genusname >> $outdir/taxa_strain.txt


