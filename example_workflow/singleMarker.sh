#!/bin/sh

#$ -N SMARKERA
#$ -hold_jid MAKECCA
#$ -t 1-250
#$ -q krt,bio,pub64

module load R

cd $1

#R --no-save --slave --args simindex.txt ccindex.txt popfile.bin $SGE_TASK_ID ccfile.bin 3000 3000 20000 logit_pvalues.txt logit_pvalues_index.txt ~/src/disease_sims/atomic_locker/atomic_locker < ~/src/disease_sims/single_marker_test.R
R --no-save --slave --args simindex.txt ccindex.txt effects.bin $SGE_TASK_ID ccfile.bin 3000 3000 20000 logit_pvalues.txt logit_pvalues_index.txt ~/src/disease_sims/atomic_locker/atomic_locker < ~/src/disease_sims/single_marker_test.R
