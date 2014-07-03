//PKG_LIBS="$LDFLAGS -lz -lboost_system" PKG_CPPFLAGS="-I$HOME/src/disease_sims $CPPFLAGS" R --no-save < CompileTest.R

#include <Rcpp.h>
#include <zlib.h>
#include <string>
#include <fstream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]		 
DataFrame getEsizes( const char * filename,
		     const unsigned long & offset )
{
  NumericVector pos,esize,count,age;
  double row[4];
  unsigned nmuts;

  //PS, zlib is fookin' magic and can read uncompressed files too...
  gzFile gzin = gzopen(filename,"rb");
  if( gzin == NULL )
    {
      Rcerr <<"Error: "
	    << filename 
	    << " could not be opened for reading\n";
      return DataFrame::create();
    }
  Rcerr << "seeking\n";
  gzseek( gzin, offset, 0 );
  gzread( gzin, &nmuts, sizeof(unsigned) );
  Rcerr << nmuts << '\n';
  for( unsigned i = 0 ; i < nmuts ; ++i )
    {
      gzread( gzin,&row[0],4*sizeof(double) );
      pos.push_back(row[0]);
      esize.push_back(row[1]);
      count.push_back(row[2]);
      age.push_back(row[3]);
    }
  gzclose(gzin);

  return DataFrame::create( Named("pos") = pos,
			    Named("esize") = esize,
			    Named("count") = count,
			    Named("age") = age );
}

// [[Rcpp::export]]	      
List getCCblock( const char * filename,
		 const unsigned long & offset )
{
  gzFile gzin = gzopen(filename,"rb");
  if( gzin == NULL )
    {
      Rcerr <<"Error: "
	    << filename 
	    << " could not be opened for reading\n";
      return List::create();
    }

  Rcerr << "seeking\n";
  gzseek( gzin,offset,0 );

  unsigned n[4];
  gzread(gzin,&n[0],4*sizeof(unsigned));

  //Read in mutation positions
  NumericVector pos(n[2]+n[3]);
  gzread(gzin,&pos[0],(n[2]+n[3])*sizeof(double));//holy crap this compiles!?!?!?

  //Matrices are allocated rows,cols
  IntegerMatrix genos(n[0]+n[1],n[2]+n[3]),burdens(n[0]+n[1],2);
  NumericMatrix phenos(n[0]+n[1],2);

  //read in the genotype matrix
  unsigned nones,ntwos;
  for( unsigned i = 0 ; i < n[0]+n[1] ; ++i )
    {
      gzread( gzin,&nones,sizeof(unsigned) );
      vector<unsigned> ones(nones);
      gzread( gzin,&ones[0],nones*sizeof(unsigned) );

      for( unsigned j = 0 ; j < nones ; ++j )
	{
	  genos(i,ones[j]) = 1;
	}

      gzread( gzin,&ntwos,sizeof(unsigned) );
      vector<unsigned> twos(ntwos);
      gzread( gzin,&twos[0],ntwos*sizeof(unsigned) );

      for( unsigned j = 0 ; j < ntwos ; ++j )
	{
	  genos(i,twos[j]) = 2;
	}
    }

  //read in burdens
  unsigned b[2];
  for( unsigned i = 0 ; i < n[0]+n[1] ; ++i )
    {
      gzread(gzin,&b[0],2*sizeof(unsigned));
      burdens(i,0)=b[0];
      burdens(i,1)=b[1];
    }

  //And phenos
  double p[2];
  for( unsigned i = 0 ; i < n[0]+n[1] ; ++i )
    {
      gzread(gzin,&p[0],2*sizeof(double));
      phenos(i,0)=p[0];
      phenos(i,1)=p[1];
    }

  gzclose(gzin);

  //Apparently, Rcpp assigns to columns first...
  //gzread( gzin,&burdens[0],2*(n[0]+n[1])*sizeof(unsigned) );
  //gzread( gzin,&phenos[0],2*(n[0]+n[1])*sizeof(double) );

  return List::create( Named("pos") = pos,
		       Named("genos") = genos,
		       Named("burdens") = burdens,
		       Named("phenos") = phenos );
}

// [[Rcpp::export]] 
NumericMatrix getPheno(const char * filename,
		       const unsigned long & offset)
{
  gzFile gzin = gzopen(filename,"rb");
  if( gzin == NULL )
    {
      Rcerr <<"Error: "
	    << filename 
	    << " could not be opened for reading\n";
      return NumericMatrix();
    }
  gzseek( gzin,offset,0 );
  unsigned n;
  gzread(gzin,&n,sizeof(unsigned));
  
  NumericMatrix rv(n,2);
  double p[2];
  for(unsigned i=0;i<n;++i)
    {
      gzread(gzin,&p[0],2*sizeof(double));
      rv(i,0)=p[0];
      rv(i,1)=p[1];
    }
  gzclose(gzin);
  return rv;
}
