genomes=$1

while read p
do
	cat "$p"/*filtered_output_final.tbl | wc -l >> filter1.txt
	cat "$p"/*filtered2_output_final.tbl | wc -l >> filter2.txt
done < $genomes