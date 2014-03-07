#ifndef __READ_CC_HPP__
#define __READ_CC_HPP__

#include <vector>
#include <limits>
struct CCblock
{

  //SN = no. segregating neutral markers
  //SC = no. segregating causative markers
  unsigned ncontrols,ncases,SN,SC;
  std::vector<double> neutral_pos,caus_pos;
  std::vector<std::vector<unsigned> > geno_matrix;
  CCblock( const unsigned & nco,
	   const unsigned & nca,
	   const std::vector<double> & np,
	   const std::vector<double> & cp,
	   const std::vector< std::vector<unsigned> > & gm );
  CCblock(void);
};

CCblock read_CC_record( const char * ccindexfile,
			const char * ccfile,
			const unsigned & recordno,
			bool * fail,
			const unsigned & max_controls = std::numeric_limits<unsigned>::max(),
			const unsigned & max_cases = std::numeric_limits<unsigned>::max(),
			const bool & rotate = false);

#endif
