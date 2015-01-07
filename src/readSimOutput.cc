#include <readSimOutput.hpp>

#include <limits>
#include <fstream>

using namespace std;

recOffsets::recOffsets(void) : pop_offset(numeric_limits<long>::max()),
			       pheno_offset(numeric_limits<long>::max()),
			       effect_offset(numeric_limits<long>::max()),
			       found(false)
{
}

recOffsets read_index( const char * idxfile,
		       const unsigned & record_no )
{
  recOffsets rv;

  ifstream index( idxfile );
  if( ! index ) { return rv; }

  unsigned ith_rep,effect_offset,pheno_offset,pop_offset;
  while( !rv.found && !index.eof() )
    {
      index >> ith_rep >> effect_offset >> pheno_offset >> pop_offset >> ws;
      if( ith_rep == record_no )
	{
	  rv.found = true;
	  rv.effect_offset=effect_offset;
	  rv.pheno_offset=pheno_offset;
	  rv.pop_offset=pop_offset;
	  break;
	}
    }
  return rv;
}

vector<effectFileData> read_effect_file( gzFile in )
{
  unsigned nrecs;
  gzread(in,&nrecs,sizeof(unsigned));

  vector< effectFileData > rv(nrecs);

  for( unsigned i = 0 ; i < nrecs ; ++i )
    {
      gzread(in,&rv[i].pos,sizeof(double));
      gzread(in,&rv[i].esize,sizeof(double));
      gzread(in,&rv[i].count,sizeof(double));
      gzread(in,&rv[i].age,sizeof(double));
    }
  return rv;
}

std::vector< std::pair<double,double> > read_phenotypes(const char * fn, z_off_t offset)
{
  std::vector< std::pair<double,double> > phenotypes;
  gzFile gzin = gzopen( fn,"rb" );
  gzseek( gzin, offset, 0);
  unsigned nphenos;
  gzread(gzin,&nphenos,sizeof(unsigned));
  for(unsigned i = 0 ; i < nphenos ; ++i )
    {
      //x is the genetic contribution to phenotype. y is the Gaussian noise from the simulation.
      //Phenotype of the individual is x+y
      double x,y;
      gzread(gzin,&x,sizeof(double));
      gzread(gzin,&y,sizeof(double)); 
      phenotypes.push_back( make_pair(x,y) );
    }
  gzclose(gzin);
  return phenotypes;
}
