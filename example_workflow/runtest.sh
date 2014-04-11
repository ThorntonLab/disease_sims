#!/bin/sh

#$ -N TESTSIMA
#$ -t 1-250
#$ -q krt,bio,pub64

module load krthornt/fwdpp/0.2.2

cd $1

SEED=`echo "$RANDOM*$SGE_TASK_ID"|bc -l`
~/src/disease_sims/TFL2013_ind -R $SGE_TASK_ID -e $2 -i simindex.txt -p popfile.bin -P phenotypes.bin -E effects.bin -S $SEED