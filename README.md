#Introduction

This repository contains various forward simulations of disease models.  More precisely, it contains forward simulations of quantitative traits under various forms of selection.

#Dependencies

1. [fwdpp](http://github.com/molpopgen/fwdpp). Version 0.2.5 or greater.
2. [libsequence](http://github.com/molpopgen/libsequence).  Version 1.7.8 or greater.
3. [boost](http://www.boost.org)
4. [GSL](http://www.gnu.org/software/gsl)

#Installation of the simulation

The project uses a configure script for setup, compiling and installation.  The complexity of the procedure depends on the complexity of your system.

##Vanilla installation instructions

For a system where you have sudo privileges and dependencie are "where they should be":

```
./configure
make
sudo make install
```
#Example workflow on UCI HPC

See the scripts in the workflow directory.

#Utilities

Used in writing the R01 for June 5, 2014

##phenoburden

Calculates statistics relatated to the deleterious mutation load in the entire population.

##CCburden

Calculates statistics related to the deleterous mutation load in a case/control panel.