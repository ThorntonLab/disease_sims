#!/bin/sh

#$ -N MAKECCA
#$ -hold_jid TESTSIMA
#$ -q krt,bio,pub64
#$ -t 1-250

module load krthornt/fwdpp/0.2.2

cd $1

SEED=$SGE_TASK_ID
~/src/disease_sims/make_case_control -i simindex.txt -p popfile.bin -P phenotypes.bin -c ccfile.bin -I ccindex.txt -r $SGE_TASK_ID -S $SEED -n 3000 -N 3000 -t 0.15 
