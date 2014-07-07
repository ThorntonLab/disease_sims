#!sh

find . -name "*.o" | xargs rm -f

echo $LDFLAGS
R_LIBS=~/R_libs_dev PKG_LIBS=$LDFLAGS PKG_CPPFLAGS="-I$HOME/src/disease_sims $CPPFLAGS" R CMD INSTALL diseaseSims
