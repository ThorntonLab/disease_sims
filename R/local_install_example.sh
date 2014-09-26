#!sh

#find . -name "*.o" | xargs rm -f

#Load the dependencies (this takes care of boost, zlib, libsequence, too, athough libsequence is not needed)
module load R
module load krthornt/fwdpp/0.2.5

#make a personal directory for R libraries

if [ ! -d ~/R_libs ]
then
    mkdir ~/R_libs
fi


R_LIBS="~/R_libs:$R_LIBS" PKG_LIBS=$LDFLAGS PKG_CPPFLAGS=$CPPFLAGS R CMD INSTALL diseaseSims
