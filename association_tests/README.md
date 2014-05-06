#Utilities for running association tests on simulation output

Utilities are a mix of custom C++ and R code.

###Installing SKAT
To install SKAT into your user's $HOME:
```
wget http://cran.r-project.org/src/contrib/SKAT_0.95.tar.gz
#You must module load first, o/w the module load will over-write what you do to the R_LIBS variable below!
module load R
tar xzf SKAT_0.95.tar.gz
cd SKAT
#Replace R_libs with whatever you are already using, if applicable:
mkdir $HOME/R_libs
export R_LIBS="$HOME/R_libs:$R_LIBS"
R CMD INSTALL -I$R_LIBS SKAT_0.95.tar.gz
```

###Loading your local installation of SKAT

```
module load R
export R_LIBS="$HOME/R_libs:$R_libs"
R
library(SKAT)
```