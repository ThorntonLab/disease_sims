#ifndef __READ_SIM_OUTPUT_HPP__
#define __READ_SIM_OUTPUT_HPP__

#include <zlib.h>
#include <vector>
#include <iostream>

struct recOffsets
{
  long pop_offset,pheno_offset,effect_offset;
  bool found;
  recOffsets(void);
};

struct effectFileData
{
  //age is in generations9
  double pos,esize,count,age;
};

recOffsets read_index( const char * idxfile,
		       const unsigned & record_no );

std::vector<effectFileData> read_effect_file( gzFile in );
#endif
