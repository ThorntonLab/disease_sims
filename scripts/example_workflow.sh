#!/usr/bin/env bash

##Cleanup from old runs
for i in simindex.txt pop.bin.gz phenotypes.bin.gz effects.bin.gz ccfile.bin.gz ccindex.txt smarkerAdditive.txt.gz smarkerRecessive.txt.gz smarkerDominant.txt.gz vexpl.pdf
do
    if [ -e $i ]
    then
	rm -f $i
    fi
done

##Run a simulation
SEED=101
TFL2013_ind -i simindex.txt -p pop.bin.gz -P phenotypes.bin.gz -E effects.bin.gz  -S $SEED -g 1000

##Make a case/control panel
make_case_control  -i simindex.txt -p pop.bin.gz -P phenotypes.bin.gz -c ccfile.bin.gz -I ccindex.txt -n 3000 -N 3000 -S $SEED -t 0.15

##Perform single-marker tests of association (logistic regression under various genetic models)
R --no-save --slave < logit.R

##Get the cumulative proportion of V_G explained due to additive effects of risk markers
R --no-save --slave < vexpl.R

##SKAT/SKAT-O
R --no-save --slave < skat.R

##burden tests (not SKAT)
