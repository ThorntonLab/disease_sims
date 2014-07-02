#include <simindex.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

typedef map<unsigned,long> mul;

//private
bool simindex::mono_increasing( const std::vector<long> & v ) const
{
  bool rv = true;
  for( std::vector<long>::const_iterator i = v.begin() + 1 ;
       rv == true && i < v.end() ; ++i )
    {
      if( !( *i > *(i-1) ) )
	{
	  rv=false;
	}
    }
  return rv;
}

//public functions
simindex::simindex(const char * filename) : e(mul()),
					    p(mul()),
					    h(mul()),
					    fileproblem(false)
{
  ifstream in(filename);
  if(!in) 
    {
      fileproblem = true;
      return;
    }
  unsigned rec;
  long ei,pi,hi;

  vector<unsigned> recs;
  vector<long> eis,pis,his;
  while(!in.eof() && !fileproblem)
    {
      in >> rec >> ei >> pi >> hi >> ws;
      if( find(recs.begin(),
	       recs.end(),
	       rec) != recs.end() ) //than the record ID exists multiple times = BAD!
	{
	  fileproblem = true;
	}

      recs.push_back(rec);
      eis.push_back(ei);
      pis.push_back(pi);
      his.push_back(hi);
    }
  in.close();

  if ( !mono_increasing(eis) || 
       !mono_increasing(pis) ||
       !mono_increasing(his) )
    {
      //index file is probably related to .gz output
      for( unsigned i = 0 ; i < recs.size() ; ++i )
	{
	  if(i == 0)
	    {
	      e[recs[i]]=eis[i];
	      p[recs[i]]=pis[i];
	      h[recs[i]]=his[i];
	    }
	  else
	    {
	      e[recs[i]]=accumulate(eis.begin(),eis.begin()+i,0l);
	      p[recs[i]]=accumulate(pis.begin(),pis.begin()+i,0l);
	      h[recs[i]]=accumulate(his.begin(),his.begin()+i,0l);
	    }
	}
    }
  else
    {
      for( unsigned i = 0 ; i < recs.size() ; ++i )
	{
	  e[recs[i]]=eis[i];
	  p[recs[i]]=pis[i];
	  h[recs[i]]=his[i];
	}
    }
}

bool simindex::file_problem() const
{
  return fileproblem;
}

bool simindex::eexists(const unsigned & i) const
{
  return ( e.find(i) != e.end() );
}

bool simindex::pexists(const unsigned & i) const
{
  return ( p.find(i) != p.end() );
}

bool simindex::hexists(const unsigned & i) const
{
  return ( h.find(i) != h.end() );
}

long simindex::eoffset(const unsigned & i) const
{
  return e.find(i)->second;
}

long simindex::poffset(const unsigned & i) const
{
  return p.find(i)->second;
}

long simindex::hoffset(const unsigned & i) const
{
  return h.find(i)->second;
}
