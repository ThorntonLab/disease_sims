#include <locking_routines.hpp>


struct flock get_whole_flock()
{
  struct flock lock;
  lock.l_type = F_WRLCK;/*Write lock*/
  lock.l_whence = SEEK_SET;
  lock.l_start = 0;
  lock.l_len = 0;/*Lock whole file*/

  return lock;
}

int open_file(const char * path,
	      int flags,
	      mode_t mode)
{
  return open(path,flags,mode);
}
