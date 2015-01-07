#!sh

rm -f simindex.txt pvalues.txt ccindex.txt esm.txt pvalues_index.txt *.bin.gz

#SEED=31588
SEED=$RANDOM
../src/TFL2013_ind -1 1000 -g 10000 -i simindex.txt -p popfile.bin.gz -P phenotypes.bin.gz -E effectsizes.bin.gz -n 0.025 -c 0.0025 -r 0.025 -e 0.1 --noise 0.075 -S $RANDOM
echo $SEED > seedfile
../src/make_case_control -i simindex.txt -p popfile.bin.gz -P phenotypes.bin.gz -c ccfile.bin.gz -I ccindex.txt -n 100 -N 100 -t 0.1

R --no-save --slave --args simindex.txt ccindex.txt effectsizes.bin.gz 0 ccfile.bin.gz 100 100 1000 pvalues.txt pvalues_index.txt  < ../R/single_marker_test.R

#./esm_chisq_zscore -c ccindex.txt -C ccfile.bin.gz -o esm.txt -r 0 -s 101 
