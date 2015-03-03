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
