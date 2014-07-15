#The diseaseSims R module

diseaseSims is an R package designed to help process simulation output, especially simulated case/control panels.

The package is a mix of R and C++ code.  The C++ code allows us to read in our simulated data from gzipped binary files and take advantage of seeking, which R's gzfile connection is unable to to reliably.

#Dependencies

The R package is not self-contained, and depends on both external libraries and a header file from the disease_sims repository.

The external dependencies are:

1. [fwdpp](https://github.com/molpopgen/fwdpp) version 0.2.4 or greater
2. [zlib](http://zlib.net), preferably version 1.2.7 or greater
3. [boost](http://boost.org), preferably 1.50 or greater
4. [R](http://www.r-project.org), 3.1 or greater
5. [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) 0.11.1 or greater

__PRO TIP:__ if you can install fwdpp, then you have dependencies 2 and 3 aready sorted out on your system, by definition.

#Installation (from source)

This is a little tricker than for a typical R package due to the dependency on a header from the simulation code.  I will assume that you are trying to install diseaseSims from the R directory of the disease_sims repository.

##Simple case: all dependent libraries are where they "should be"

If the first three dependencies listed above are in standard locations (typically somewhere in /usr or /usr/local), then the following command will suffice to install the R library:

```
PKG_CPPFLAGS=$CPPFLAGS R CMD INSTALL diseaseSims
```

##More complex case: your libraries are in funny places and/or your system uses [modules](http://modules.sourceforge.net/) to make things like libraries available to users

Basically, if you have libraries in your home directory, or various other places where the C/C++ compilers do NOT look by default, then you typically need to adjust both CPPFLAGS and LDFLAGS to tell the preprocessor and linker, respectively, where to look.

For the sake of argument, let's assume that dependencies are in your user's $HOME directory.  Specifically, headers are in $HOME/include and run-time libraries are in $HOME/lib.  Unless you like going insane, you probably have these variables already set in a file like your .bashrc:

```
CPPFLAGS=-I$HOME/include 
LDFLAGS=-L$HOME/lib
```

OK, we can work with that:

```
PKG_LIBS=$LDFLAGS PKG_CPPFLAGS=$CPPFLAGS R CMD INSTALL diseaseSims
```

Basically, PKG_LIBS is a combo of LDFLAGS plus the individual -l flags uses to link libraries (those flags are already included in the diseaseSims Makevars file, so we don't have to specify them here).

If you use the modules system, then the "module load XXXX" commands should be prepending the correct -L and -I paths to LDFLAGS and CPPFLAGS, respectively.  If so, then the above command will work for you.  If not, talk to your sysadmin because doing this right is the whole point of the module system.

##Installation into user's $HOME/R_libs

```
R_LIBS=$HOME/R_libs:$R_LIBS PKG_LIBS=$LDFLAGS PKG_CPPFLAGS=$CPPFLAGS R CMD INSTALL diseaseSims
```

###Installation on the UCI HPC into a user's $HOME

Execute the local_install_example.sh in disease_sims/R.