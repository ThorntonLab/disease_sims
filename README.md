#Introduction

This repository contains various forward simulations of disease models.  More precisely, it contains forward simulations of quantitative traits under various forms of selection.

#Dependencies

1. [fwdpp](http://github.com/molpopgen/fwdpp). Version 0.2.4 or greater required, but 0.2.5 or greather _highly_ recommended for performance reasons.
2. [libsequence](http://github.com/molpopgen/libsequence).  Version 1.7.8 or greater.
3. [boost](http://www.boost.org)
4. [GSL](http://www.gnu.org/software/gsl)

#Installation of the simulation

The project uses a configure script for setup, compiling and installation.  The complexity of the procedure depends on the complexity of your system.

##Vanilla installation instructions

For a system where you have sudo privileges and dependencies are "where they should be":

```
./configure
make
sudo make install
```

If you have dependencies in locations other than /usr/local.  For example, dependencies are in $HOME

```
./configure CXXFLAGS=-I$HOME/include LDFLAGS=-L$HOME/lib
```

To install the simulation in your user's $HOME:

```
./configure --prefix=$HOME
```

The above may be mixed and matched as needed.

##Notes

On many systems, LD_LIBRARY_PATH may be needed to run make_case_control.

#The simulation programs

##TFL2013_ind

This is a complete rewrite of the program used in our previous paper (Thornton _et al._. 2013, doi:10.1371/journal.pgen.1003258).  If you run the program without arguments, you'll get an extensive printout to screen about its usage.  This section of this document will attempt to put that information into context.

###The demographic model

The program simulates a single population of _N_ diploids for _n1_ generations.  Optionally, the simulation may continue for another _n2_ generations, during which the population size changes (exponentially) from _N_ to _N2_.

For example, to simulation a population of _N=10,000_ for 50,000 generations, and then grow it to _N=100,000_ over a period of 10,000 generations, say

~~~
-1 10000 -g 50000 -2 100000 -G 10000
~~~

Let's look at this set of options:

~~~
-1 10000 -g 50000 -2 100000
~~~

Well, that's close, but the default length of the growth period is 0, so no growth will happen and the simulation will end after 50,000 generations.

###Global parameters of a locus

There are several parameters that are used for all genetic models:

* The mutation rate to variants not affecting phenotype/fitness (hereafter "neutral" mutations arising at the neutral mutation rate), specified by --neutral/-n and a floating-point value > 0.
* The mutation rate to variants affecting phenotype/fitness (hereafter "causative" mutations arising at the deleterious mutation rate), specified by --causative/c and a floataing-point value > 0.
* The mean effect size of a causative mutation, specified by --esize/-e and a floating point value > 0.
* The recombination rate (per diploid, per generation), specified by --recrate/-r and a floating point value > 0.

Example:

~~~
-n 0.00125 -c 0.000125 -e 0.1
~~~

The mutation rates are _per gamete, per generation_, and are taken as the mean of a Poisson distribution.  The scale of the value pased to --esize/-e depends on the specific genetic model (more on this below).

####The model of a locus

In this simulation, a "gene region" is modeled as a series of positions on the continuous, half-open interval (0,1].  Mutations and crossover events occur uniformly along this region.  Mutations occur according to the infinitely-many sites model, meaning that each new mutation in the population occurs at a novel position.

###Output files

The program generates up to three output files plus an index file recording the sizes of various records. 

* -i/--indexfile _filename_ specifies the name of the index file.
* --popfile/-p _filename_ specifies an output file where the entire population will be written.
* --effectsfile/-E _filename_ specifies and output file where some information about segregating mutations will be written.
* --phenotypes/-O _filename_ specifies an output file where individual diploid phenotypes will be written. (_This file is not written for all genetic models._)

Only the first three are mandatory.

#### The index file

The index file contains rows of 4 integer columns:

1. The replicate ID number, specified by --replicate/-R and an integer >= 0
2. The size of the data written to the "effects file", in bytes.
3. The size of the data written to the "phenotypes file", in bytes.  (A value of -1 is used when the phenotypes file is not written.)
4. The size of the data written to the "population file", in bytes.

The index file is plain-text, and human-readable.  The remaining files are not human-readable.  Their data are written in a native binary format and the files them selves are gzip-compressed.

#### The phenotypes file

For each simulated replicate, the phenotypes file contains the following information:

1. A 32-bit unsigned integer representing the number of diploids present in the replicate.  Let's call this number _D_.
2. _2D_ floating point values whose size is equal to sizeof(double) on the machine where the simulation was compiled.  Each pair of doubles is the genetic (_G_) and random (_E_) component of phenotype.

A record may be read into R using the readBin function:

~~~{r}
#open for reading in binary mode
f=file("phenotypes.bin.gz","rb");
D=readBin(f,"integer",1)
m=matrix(readBin(f,"numeric",2*D),ncol=2,byrow=TRUE)
~~~

The heritability due to the locus for this replicate is simply:

~~~{r}
h2 = var( m[,1] )/( var(m[,1]) + var(m[,2]) )
~~~

You may read all the replicates in from a file using a for loop, etc.  Please note that seeking within a gzipped file doesn't work all that well in R.


#Example workflow on UCI HPC

See the scripts in the workflow directory.

#Utilities

Used in writing the R01 for June 5, 2014

##phenoburden

Calculates statistics relatated to the deleterious mutation load in the entire population.

##CCburden

Calculates statistics related to the deleterous mutation load in a case/control panel.


#Important links:

* PCGC software: https://sites.google.com/site/davidgolanshomepage/software/pcgc
