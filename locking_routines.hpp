#ifndef __LOCKING_ROUTINES_HPP__
#define __LOCKING_ROUTINES_HPP__

#include <fcntl.h>

struct flock get_whole_flock(void);

int open_file(const char * path,
	      int flags = O_APPEND|O_CREAT|O_WRONLY,
	      mode_t mode = S_IRWXU | S_IRWXG);

#endif
