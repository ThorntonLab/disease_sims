#Introduction

This repository contains various forward simulations of disease models.  More precisely, it contains forward simulations of quantitative traits under various forms of selection.

#Dependencies

1. [fwdpp](http://github.com/molpopgen/fwdpp). Version 0.2.5 or greater.
2. [libsequence](http://github.com/molpopgen/libsequence).  Version 1.8.4 or greater.
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
