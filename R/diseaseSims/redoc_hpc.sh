#!sh

module load R krthornt/fwdpp/0.2.4
R_LIBS=~/$HOME/R_libs:$R_LIBS PKG_LIBS=$LDFLAGS PKG_CPPFLAGS="-I$HOME/src/disease_sims $CPPFLAGS" R_LIBS=/data/users/krthornt/R_libs:$R_LIBS sh redoc.sh