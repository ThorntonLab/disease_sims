#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <zlib.h>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

using namespace Rcpp;
using namespace std;

struct boostScopedLockManager
{
  typedef boost::interprocess::file_lock flock_t;
  typedef boost::interprocess::scoped_lock<flock_t> slock_t;
  gzFile gzf;
  flock_t * flock;
  slock_t * slock;
  bool closed;
  
  void finish()
  {
    if(closed) return;
    flock->unlock();
    gzclose(gzf);
    closed = true;
  }
  
  boostScopedLockManager( const char * filename,
			  const char * mode) :
    gzf(gzopen(filename,mode)),
    flock(new flock_t(filename)),
    slock(new slock_t(*flock)),
    closed(false)
  {
  }

  ~boostScopedLockManager()
  {
    flock->unlock();
    slock->unlock();
    if(!closed)
      {
	gzclose(gzf);
      }
    delete flock;
    delete slock;
  }
};

//' Manage access to a gz-compressed file
//' @param filename The file name to write to
//' @param mode The mode for opening the file.
//' @return An external pointer to a boostScopedLockManager
//[[Rcpp::export("gzLock")]]
SEXP initZlibBoostScopedLock( const char * filename,
			      const char * mode )
{
  return XPtr<boostScopedLockManager> ( new boostScopedLockManager(filename,mode) );
}

//' Unlock and close a gz-compressed file
//' @param s A boostScopedLockManager
//[[Rcpp::export("gzUnlock")]]
void endZlibBoostScopedLock( SEXP s )
{
  XPtr<boostScopedLockManager> ps(s);
  ps->finish();
}

//' Write a data frame to file
//' @param s A data frame
//' @param locker A boostScopedLockManager
//' @param colnames If TRUE, write column names to file
//' @param sep column separater
//[[Rcpp::export]]
int writeDataFrame( const SEXP s,
		    SEXP locker,
		    const bool & colnames = true,
		    const char * sep = "\t")
{
  List df(s);
  XPtr<boostScopedLockManager> lfile(locker);
  int ncol = df.size();
  if(!ncol) return 0;
  int nrow = CharacterVector(df[0]).length();
  CharacterVector names = df.attr("names");

  ostringstream o;

  if(colnames)
    {
      int c = 0;
      for( ; c < ncol-1 ; ++c)
  	{
  	  o << names[c] << sep;
  	}
      o << names[c] << '\n';
    }
  for(int r = 0; r < nrow ; ++r )
    {
      int c = 0;
      for( ; c < ncol-1 ; ++c )
  	{
  	  o << CharacterVector(df[c])[r] << sep;
  	}
      o << CharacterVector(df[c])[r] << '\n';
    }

  return gzwrite(lfile->gzf,o.str().c_str(),o.str().size());
}

//' Write a matrix to file
//' @param s A matrix
//' @param locker A boostScopedLockManager
//' @param colnames If TRUE, write column names to file
//' @param sep column separater
//[[Rcpp::export]]
int writeMatrix( const SEXP s,
		 SEXP locker,
		 const bool & colnames = true,
		 const char * sep = "\t")
{
  return writeDataFrame(DataFrame(s),locker,colnames,sep);
}
