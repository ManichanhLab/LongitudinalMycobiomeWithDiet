SparCC.py merged_func.csv -c cor.txt -v cov.txt
MakeBootstraps.py merged_func.csv -p simulated/ -t permuted_#
mkdir boot_cor
mkdir boot_cov
for i in `seq 0 99`; do SparCC.py simulated/permuted_$i -c boot_cor/simulated_$i.txt -v boot_cov/simulated_$i.txt >> log; done
mkdir pvals
PseudoPvals.py cor.txt boot_cor/simulated_#.txt 100 -o pvals/one_sided.txt -t one_sided
