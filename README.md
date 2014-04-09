#Introduction

This repository contains various forward simulations of disease models.  More precisely, it contains forward simulations of quantitative traits under various forms of selection.

#Dependencies

1. [fwdpp](http://github.com/molpopgen/fwdpp). Version 0.2.0 or greater.
2. [libsequence](http://github.com/molpopgen/libsequence).  Version 1.7.8 or greater.
3. [boost](http://www.boost.org)
4. [GSL](http://www.gnu.org/software/gsl)

#Example workflow on UCI HPC

##Script to run the simulation

```
#!/bin/sh

#$ -N TESTSIM
#$ -t 1-250
#$ -q krt,bio,pub64

module load krthornt/fwdpp/0.2.2

cd /fast-scratch/krthornt/test_disease_sim

SEED=`echo "$RANDOM*$SGE_TASK_ID"|bc -l`
~/src/disease_sims/TFL2013_ind -R $SGE_TASK_ID -e 0.075 -i simindex.txt -p popfile.bin -P phenotypes.bin -E effects.bin -S $SEED
```

##Script to make case/control panels
```
#!/bin/sh

#$ -N MAKECC
#$ -hold_jid TESTSIM
#$ -q krt,bio,pub64
#$ -t 1-250

cd /fast-scratch/krthornt/test_disease_sim

SEED=`echo "$RANDOM*$SGE_TASK_ID"|bc -l`
~/src/disease_sims/make_case_control -i simindex.txt -p popfile.bin -P phenotypes.bin -c ccfile.bin -I ccindex.txt -r $SGE_TASK_ID -S $SEED -n 3000 -N 3000 -t 0.15
```

##Perform single-marker logistic regression of case/control status onto genotype under an additive model
```
#!/bin/sh

#$ -N SMARKER
#$ -hold_jid MAKECC
#$ -t 1-250

module load R

cd /fast-scratch/krthornt/test_disease_sim

R --no-save --slave --args simindex.txt ccindex.txt popfile.bin $SGE_TASK_ID ccfile.bin 3000 3000 20000 logit_pvalues.txt logit_pvalues_index.txt ~/src/disease_sims/atomic_locker/atomic_locker < ~/src/disease_sims/single_marker_test.R
```

##Perform the ESM test over a variety of "K" values
```
#!/bin/sh

#$ -N ESM
#$ -hold_jid MAKECC
#$ -t 1-250

module load krthornt/fwdpp/0.2.2

cd /fast-scratch/krthornt/test_disease_sim

SEED=`echo "$RANDOM*$SGE_TASK_ID"|bc -l`

for K in 50 150 250
do
    ~/src/disease_sims/esm_chisq_zscore -c ccindex.txt -C ccfile.bin -o esm.K$K.txt -r $SGE_TASK_ID -s $SEED -K $K
done
```
