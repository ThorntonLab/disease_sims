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

##Installing the R package

Once the main simulation is installed, you may install the diseaseSims R package.  The R package requires headers and a shared library installed when the simulation is installed.

The R package has the following additional dependencies:

1. [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html)
2. [dplyr](http://cran.r-project.org/web/packages/dplyr/index.html)
3. [buRden](http://github.com/molpopgen/buRden)

The following packages are suggested:

1. [SparseM](http://cran.r-project.org/web/packages/SparseM/index.html)

###Vanilla instructions for installing the R package

Once the above dependencies are accounted for,

```
cd R
R CMD INSTALL diseaseSims
```

### Build the documentation for the R package

```
R CMD Rd2pdf diseaseSims
```

##Notes

On many systems, LD_LIBRARY_PATH may be needed to run make_case_control.

# Quick-start guide for the impatient

This section assumes a bash/bash-like shell.  Further, the command lines shown here are explicit, and many of the options specified below are actually just the default values.

Simulate a population of N = 2,000 diploids for 8N generations, and then grow to N = 100,000 in 500 generation.  The genetic model and locus parameters are taken from Thornton et al. (2013), and the mean effect size of a risk mutation is 0.1:

~~~
TFL2013_ind -1 2000 -g 16000 -2 100000 -G 500 -n 0.00125 -c 0.000125 -r 0.00125 -e 0.1 --noise 0.075 -i simindex.txt -p popfile.bin.gz -P phenotypes.bin.gz -E effects.bin.gz -S $RANDOM
~~~

Make a case/control panel of 3,000 controls + 3,000 cases.  Cases are from the upper 10% of the distribution of phenotypes in the population, and controls are sampled from the mean +/- 0.5 standard deviations.

~~~
make_case_control -i simindex.txt -p popfile.bin.gz -P phenotypes.bin.gz -I ccindex.txt -c ccfile.bin.gz -n 3000 -N 3000 -t 0.10 --control-range 0.5 --seed $RANDOM
~~~

A lot of the other types of post-simulation processing are handled by the R package.  See the vignettes, manual, etc., for help.

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

###The model of a locus

In this simulation, a "gene region" is modeled as a series of positions on the continuous, half-open interval (0,1].  Mutations and crossover events occur uniformly along this region.  Mutations occur according to the infinitely-many sites model, meaning that each new mutation in the population occurs at a novel position.

###Output files

The program generates up to three output files plus an index file recording the sizes of various records. 

* -i/--indexfile _filename_ specifies the name of the index file.
* --popfile/-p _filename_ specifies an output file where the entire population will be written.
* --effectsfile/-E _filename_ specifies and output file where some information about segregating mutations will be written.
* --phenotypes/-P _filename_ specifies an output file where individual diploid phenotypes will be written. (_This file is not written for all genetic models._)

Only the first three are mandatory.

#### The index file

The index file contains rows of 4 integer columns:

1. The replicate ID number, specified by --replicate/-R and an integer >= 0
2. The size of the data written to the "effects file", in bytes.
3. The size of the data written to the "phenotypes file", in bytes.  (A value of -1 is used when the phenotypes file is not written.)
4. The size of the data written to the "population file", in bytes.

The index file is plain-text, and human-readable.  The remaining files are not human-readable.  Their data are written in a native binary format and the files themselves are gzip-compressed.

#### The phenotypes file

For each simulated replicate, the phenotypes file contains the following information:

1. A 32-bit unsigned integer representing the number of diploids present in the replicate.  Let's call this number _D_.
2. _2D_ floating point values whose size is equal to sizeof(double) on the machine where the simulation was compiled.  Each pair of doubles is the genetic (_G_) and random (_E_) component of phenotype.

A record may be read into R using the readBin function:

~~~{r}
#open for reading in binary mode
f=gzfile("phenotypes.bin.gz","rb");
D=readBin(f,"integer",1)
m=matrix(readBin(f,"numeric",2*D),ncol=2,byrow=TRUE)
~~~

The heritability due to the locus for this replicate is simply:

~~~{r}
h2 = var( m[,1] )/( var(m[,1]) + var(m[,2]) )
~~~

You may read all the replicates in from a file using a for loop, etc.  Please note that seeking within a gzipped file doesn't work all that well in R.

#### The effect sizes file

For each simulated replicate, the following data are recorded:

1. A 32-bit unsigned integer representing the number of mutations segregating in the replicate. Let's call this number _M_.
2. There are then _M_ sets of 4 doubles.  These doubles are the mutation position, its effect size (or selection coefficient), its number of occurrences in the population, and its age in generations.

Again, this file is mean for easy processing in R:

~~~{r}
f=gzfile("effectsfile.bin.gz","rb")
M=readBin(f,"integer",1)
muts = matrix(readBin(f,"numeric",4*M),ncol=4,byrow=TRUE)
~~~

#### The population file.

This file contains everything about the currently-segregating variation in the population.  Its format is complex, and not readable in R, python, etc., because it contains "deep copies" of internal [fwdpp](http://molpopgen.github.io/fwdpp/) data structures.  If you want to write programs accessing the data in this file, please see the following resources:

1. The fwdpp reference manual, found [here](http://molpopgen.github.io/fwdpp/).  Take a look at "Related Pages", and then "Tutorial 3: Data serialization".
2. The first several lines of the main routine of make_case_control.cc, in the src subdirectory of this repository.  People with C++ experience should look for this function call:

~~~{cpp}
read_binary_pop( &gametes, &mutations, &diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
~~~

When you understand what goes into that, you'll much of what you need to be able to process the contents of this file.

### File locking: you do not need to write each replicate to a separate file!

The program implements something called "file locking", which allows multiple independent processes to share the same output file.  If you want to run a bunch of replicates with the same parameters but different seeds, you simply need to change the replicated ID number (--replicate/-R) and specify the same output file for each process.  When a replicate is done, the program will talk to the system and request access to the output file.  If it gets access, it "locks" the file, and appends its data to the end.  If another process has locked the file, it'll wait until that lock is released.

This file locking makes the output of the simulation much more manageable, and it is very easily used on cluster systems supporting "array jobs".

### Genetic models implemented

The program implements a suite of different models of gene action ("genetic models").

#### "TFL2013"

This is the model of non-complementing recessive mutations that we introduced in Thornton _et al._. 2013, doi:10.1371/journal.pgen.1003258.  This is the default genetic model, and you can see the quantitative details of the model in the PLoS Genetics paper.

The relevant program parameters and their intrepretation are:

Parameter [value] | Interpretation
--------- | ------------
-e/--esize [positive double] | The effect size of a causative mutation.  By default, this is the mean of an exponential distribution ("lambda" in the PLoS Genetics paper).  If -C/--constant is used, the effect sizes are fixed at the input value
--noise [positive double] | The standard deviation of Gaussian noise added to phenotype (the mean of the Gaussian equals 0).  The default is the value used in the TFL2013 paper.  This parameter is used to tune the heritability.  A value of 0 will make the heritability equal to 1.  The trait value of a diploid is P = G + E, where G is the geometric mean of maternal and paternal haplotype effect sizes, and E is a Gaussian noise term with mean zero.
--sigma [positive double] | The standard deviation of the Gaussian fitness function.  Leave this at 1 (the default) unless you want to go insane.
--optimum [double] | For growth models, change the optimum trait value from 0 to the user-input value.  This changes the Gaussian fitness function from N(0,sigma) to N(optimum,sigma)

#### Additive and multiplicative models.

I lump these two models together as the former approximates the latter very closely for small effect sizes.

Model | Description
------ | ----------
--additive | Change the genetic model to additive across causative mutations.  For this model, G becomes the sum of the effect sizes of all causative mutations in a diploid.
--multiplicative | Change the genetic model to multiplicative across causative mutations.  For this model, G becomes the (product of (1 + e_i)) - 1, where e_i is the effect size of the i-th mutation.  This is Risch's model, with an additional - 1 subtracted from G such that the trait value in the absence of causative mutations equals 0.


The parameters for this model are the same as the above:

Parameter [value] | Interpretation
--------- | ------------
-e/--esize [positive double] | The effect size of a causative mutation.  By default, this is the mean of an exponential distribution ("lambda" in the PLoS Genetics paper).  If -C/--constant is used, the effect sizes are fixed at the input value
--noise [positive double] | The standard deviation of Gaussian noise added to phenotype (the mean of the Gaussian equals 0).  The default is the value used in the TFL2013 paper.  This parameter is used to tune the heritability.  A value of 0 will make the heritability equal to 1.  The trait value of a diploid is P = G + E, where G is the geometric mean of maternal and paternal haplotype effect sizes, and E is a Gaussian noise term with mean zero.
--sigma [positive double] | The standard deviation of the Gaussian fitness function.  Leave this at 1 (the default) unless you want to go insane.
--optimum [double] | For growth models, change the optimum trait value from 0 to the user-input value.  This changes the Gaussian fitness function from N(0,sigma) to N(optimum,sigma)

#### Multiplicative with dominance

For the multiplicative model above, G = (product of (1 + e_i)) - 1.  However, it is common in population genetics to allow the dominance of a causative mutation to vary.  If we consider the following 

Genotype | Meaning | Contribution to trait value
-------  | ------- | ----------
AA | Homozygous non-mutant | 1
Aa | Heterozygous mutant | 1 + h*e
aa | Homozygous mutant | 1 + 2*e

For the model with dominance, G = (product of (1+g_i)) -1, where g_i is the genotype at the i-th position containing a causative mutation.

Parameter [value] | Interpretation
--------- | ------------
-e/--esize [positive double] | The effect size of a causative mutation.  By default, this is the mean of an exponential distribution ("lambda" in the PLoS Genetics paper).  If -C/--constant is used, the effect sizes are fixed at the input value
-d/--dominance [positive double] | The dominance of a causative mutation
--noise [positive double] | The standard deviation of Gaussian noise added to phenotype (the mean of the Gaussian equals 0).  The default is the value used in the TFL2013 paper.  This parameter is used to tune the heritability.  A value of 0 will make the heritability equal to 1.  The trait value of a diploid is P = G + E, where G is the geometric mean of maternal and paternal haplotype effect sizes, and E is a Gaussian noise term with mean zero.
--sigma [positive double] | The standard deviation of the Gaussian fitness function.  Leave this at 1 (the default) unless you want to go insane.
--optimum [double] | For growth models, change the optimum trait value from 0 to the user-input value.  This changes the Gaussian fitness function from N(0,sigma) to N(optimum,sigma)

#### Eyre-Walker's model

Model | Description
------ | -------
--eyre-walker | Change the genetic model to that of Eyre-Walker (2010), www.pnas.org/cgi/doi/10.1073/pnas.0906182107.

For this model, no phenotypes are actually simulated.  Rather, the simulation is performed with selection under an additive model.  Thus, the usage of the program changes somewhat.

For this model, scaled selection coefficients (4Ns) are drawn from a Gamma distribution with shape parameter __Beta__ and mean __S__.

Parameter [value] | Interpretation
--------- | ------------
-e/--esize [positive double] | This value will be used for __S__.
--ewshape [positive double] | This value will be used for __Beta__.
--noise [positive double] | Has no effect.
--sigma [positive double] | Has no effect.
--optimum [double] | Has no effect
--phenotypes/-P [string] | Has no effect.  No trait values are actually simulated, so there are no phenotypes to write out.

What does the scaled selection parameter 4Ns mean in a growing population?  If we consider the demographic model described above, we have a population growing from size _N1_ to size _N2_ in _G_ generations.  Given that the simulation works forwards in time, we have implemented the most "natural" definition of 4Ns, which is to take it to refer to 4N1s.  In other words, selection coefficients are scaled __with resepct to the ancestral population size__ in growth models.  Be very careful when comparing results of different simulations, as some may take 4Ns to mean 4N2s!!!!!!

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
