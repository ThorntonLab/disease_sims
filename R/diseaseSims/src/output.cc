#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

using namespace Rcpp;
using namespace std;

struct boostScopedLockManager
{
  typedef boost::interprocess::file_lock flock_t;
  typedef boost::interprocess::scoped_lock<flock_t> slock_t;
  int existed,fd;
  flock_t flock;
  slock_t slock;
  bool closed,unlocked;
  
  void finish()
  {
    if( closed ) return;
    flock.unlock();
    slock.unlock();
    unlocked = true;
    if( existed ) close(fd);
    closed = true;
  }
  
  boostScopedLockManager( const char * filename ) :
    existed(1),fd(0),flock(flock_t()),slock(slock_t()),closed(false),unlocked(false)
  {
    //1. Check if the file exists
    struct stat buffer;
    existed = stat(filename,&buffer);

    if( existed != 0 ) //then the file does not exists
      {
	//Create the file
	fd = open(filename, O_WRONLY|O_CREAT|O_APPEND,S_IRUSR|S_IWUSR);
	if( fd == -1 )
	  {
	    ostringstream o;
	    o << "could not create file descriptor for " << filename
	      << " with flags O_WRONLY|O_CREAT|O_APPEND" 
	      << " and permissions S_IRUSR|S_IWUSR (user read/write)."
	      << " Line " << __LINE__ << " of " << __FILE__ << '\n';
	    Rcpp::stop(o.str());
	  }
      }
    flock = flock_t(filename);
    slock = slock_t(flock);
  }

  ~boostScopedLockManager()
  {
    if(!unlocked)
      {
	flock.unlock();
	slock.unlock();
      }
    if(!closed&&existed)
      {
	close(fd);
      }
  }
};

//' Manage access to a file
//' @param filename The file name to write to
//' @return An external pointer to a boostScopedLockManager
//[[Rcpp::export("Lock")]]
SEXP initScopedLock( const char * filename )
{
  return XPtr<boostScopedLockManager> ( new boostScopedLockManager(filename) );
}

//' Unlock and close a file
//' @param s A boostScopedLockManager
//' @note Failing to call this function may cause file locks to persist until the R session ends
//[[Rcpp::export("Unlock")]]
void endScopedLock( SEXP s )
{
  XPtr<boostScopedLockManager> ps(s);
  ps->finish();
}

