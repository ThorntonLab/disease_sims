#!/bin/sh

#$ -q krt,bio,pub64
#$ -t 1-250
module load R

export R_LIBS="$HOME/R_libs:$R_LIBS"

/usr/bin/time -f "%e %M" -o SKATtest_time2.txt R --no-save --slave --args ccindex.txt ccfile.bin $SGE_TASK_ID foo_index.txt foo 5000 5000 gwas ~/src/disease_sims/atomic_locker/atomic_locker < ~/src/disease_sims/association_tests/runSKAT.R

