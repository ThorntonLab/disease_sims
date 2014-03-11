#!sh

rm -f simindex.txt pvalues.txt ccindex.txt pvalues_index.txt *.bin

./TFL2013_ind -1 1000 -g 10000 -i simindex.txt -p popfile.bin -P phenotypes.bin -E effectsizes.bin -S 101
./make_case_control -i simindex.txt -p popfile.bin -P phenotypes.bin -c ccfile.bin -I ccindex.txt -n 100 -N 100 -t 0.1

R --no-save --args simindex.txt ccindex.txt effectsizes.bin 0 ccfile.bin 100 100 1000 pvalues.txt pvalues_index.txt ./atomic_locker/atomic_locker < single_marker_test.R
