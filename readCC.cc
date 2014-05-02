#include <readCC.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

CCblock::CCblock( const unsigned & nco,
		  const unsigned & nca,
		  const vector<double> & np,
		  const vector<double> & cp,
		  const vector< vector<unsigned> > & gm,
		  const vector< std::pair<unsigned,unsigned> > & b,
		  const vector< std::pair<double,double> > & p ) : ncontrols(nco), ncases(nca), 
								   SN(np.size()),SC(cp.size()),
								   neutral_pos(np), caus_pos(cp),
								   geno_matrix(gm),
								   burden(b),
								   phenos(p)
{
}

CCblock::CCblock(void) : ncontrols(0), ncases(0),
			 neutral_pos( vector<double>() ),
			 caus_pos( vector<double>() ),
			 geno_matrix( vector<vector<unsigned> >() )
{
}

CCblock read_CC_record( const char * ccindexfile,
			const char * ccfile,
			const unsigned & recordno,
			bool * fail,
			const unsigned & max_controls,
			const unsigned & max_cases,
			const bool & rotate)
{
  ifstream in(ccindexfile);
  unsigned r;
  unsigned long long pos;
  *fail = true;
  while(!in.eof() && *fail)
    {
      in >> r >> pos >> ws;
      if(r==recordno)
	{
	  *fail=false;
	}
    }
  in.close();

  in.open(ccfile,ios_base::in | ios_base::binary);
  in.seekg(pos);

  unsigned ncont,ncase,SN,SC;

  in.read( reinterpret_cast<char *>(&ncont),sizeof(unsigned) );
  in.read( reinterpret_cast<char *>(&ncase),sizeof(unsigned) );
  in.read( reinterpret_cast<char *>(&SN),sizeof(unsigned) );
  in.read( reinterpret_cast<char *>(&SC),sizeof(unsigned) );

  vector<double> npos(SN),cpos(SC);
  vector<vector<unsigned> > genos( (rotate == false) ? min(ncont,max_controls)+min(ncase,max_cases) : (SN + SC),
				   (rotate == false) ? vector<unsigned>(SN + SC) :  vector<unsigned>(min(ncont,max_controls)+min(ncase,max_cases))
				   );
  in.read( reinterpret_cast<char *>(&npos[0]), SN*sizeof(double) );
  in.read( reinterpret_cast<char *>(&cpos[0]), SC*sizeof(double) );

  unsigned NCONT=min(ncont,max_controls),NCASE=min(ncase,max_cases);
  unsigned conts_stored=0,cases_stored=0,nm;
  for( unsigned i=0; i<ncont+ncase ; ++i )
    {
      vector<unsigned> buffer(SN+SC,0u),idx;
      //Read in genotype values of 1
      in.read( reinterpret_cast<char *>(&nm),sizeof(unsigned) );
      idx.resize(nm);
      in.read( reinterpret_cast<char *>(&idx[0]),nm*sizeof(unsigned) );
      for( unsigned j = 0 ; j < nm ; ++j )
	{
	  buffer[ idx[j] ] = 1;
	}
      //followed by genotype values of 2
      in.read( reinterpret_cast<char *>(&nm),sizeof(unsigned) );
      idx.resize(nm);
      in.read( reinterpret_cast<char *>(&idx[0]),nm*sizeof(unsigned) );
      for( unsigned j = 0 ; j < nm ; ++j )
	{
	  buffer[ idx[j] ] = 2;
	}
      //in.read( reinterpret_cast<char *>(&buffer[0]),(SN+SC)*sizeof(unsigned) );
      if( i < ncont && conts_stored < NCONT) 
	{
	  if(!rotate)
	    {
	      genos[i]=buffer;
	    }
	  else
	    {
	      for(unsigned j = 0 ; j < buffer.size() ; ++j)
		{
		  genos[j][i] = buffer[j];
		}
	    }
	  ++conts_stored;
	}
      else if (i >= ncont && cases_stored < NCASE) 
	{
	  if(!rotate)
	    {
	      genos[i]=buffer;
	    }
	  else
	    {
	      for(unsigned j = 0 ; j < buffer.size() ; ++j)
		{
		  genos[j][i] = buffer[j];
		}
	    }
	  ++cases_stored;
	}
      else //skip the record
	{
	}
    }
  //read in the number of causative mutations on each haplo
  vector< std::pair<unsigned,unsigned> > burdens;
  unsigned b[2];
  for( unsigned i=0; i<ncont+ncase ; ++i )
    {
      in.read( reinterpret_cast<char *>(&b[0]),2*sizeof(unsigned) );
      burdens.push_back( std::make_pair(b[0],b[1]) );
    }
  //read in the number of causative mutations on each haplo
  vector< std::pair<double,double> > phenos;
  double p[2];
  for( unsigned i=0; i<ncont+ncase ; ++i )
    {
      in.read( reinterpret_cast<char *>(&p[0]),2*sizeof(double) );
      phenos.push_back( std::make_pair(p[0],p[1]) );
    }
  assert( conts_stored == NCONT );
  assert( cases_stored == NCASE );
  return CCblock( NCONT,NCASE,npos,cpos,genos,burdens,phenos );
}
