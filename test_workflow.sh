#!sh

rm -f simindex.txt pvalues.txt ccindex.txt esm.txt pvalues_index.txt *.bin

./TFL2013_ind -1 1000 -g 10000 -i simindex.txt -p popfile.bin -P phenotypes.bin -E effectsizes.bin -n 0.025 -c 0.0025 -r 0.025 -e 0.1 --noise 0.075 -S $RANDOM
./make_case_control -i simindex.txt -p popfile.bin -P phenotypes.bin -c ccfile.bin -I ccindex.txt -n 100 -N 100 -t 0.1

R --no-save --slave --args simindex.txt ccindex.txt effectsizes.bin 0 ccfile.bin 100 100 1000 pvalues.txt pvalues_index.txt ./atomic_locker/atomic_locker < single_marker_test.R

./esm_chisq_zscore -c ccindex.txt -C ccfile.bin -o esm.txt -r 0 -s 101 
