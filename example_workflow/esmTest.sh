#!/bin/sh

#$ -N ESMA
#$ -hold_jid MAKECCA
#$ -t 1-250
#$ -q krt,bio,pub64

module load krthornt/fwdpp/0.2.2

cd $1

SEED=`echo "$RANDOM*$SGE_TASK_ID"|bc -l`

for K in 50 150 250
do
    ~/src/disease_sims/esm_chisq_zscore -c ccindex.txt -C ccfile.bin -o esm.K$K.txt -r $SGE_TASK_ID -s $SEED -K $K
done