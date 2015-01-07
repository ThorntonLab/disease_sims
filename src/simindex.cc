#include <simindex.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace std;

typedef map<unsigned,long> mul;


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
  long ei=0,pi=0,hi=0,ei_t,pi_t,hi_t;

  vector<unsigned> recs;
  vector<long> eis,pis,his;
  while(!in.eof() && !fileproblem)
    {
      in >> rec >> ei_t >> pi_t >> hi_t >> ws;
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
      
      ei += ei_t;
      pi += pi_t;
      hi += hi_t;
    }
  in.close();

  for( unsigned i = 0 ; i < recs.size() ; ++i )
    {
      e[recs[i]]=eis[i];
      p[recs[i]]=pis[i];
      h[recs[i]]=his[i];
    }
  // if ( !mono_increasing(eis) || 
  //      !mono_increasing(pis) ||
  //      !mono_increasing(his) )
  //   {
  //     //index file is probably related to .gz output
  //     for( unsigned i = 0 ; i < recs.size() ; ++i )
  // 	{
  // 	  if(i == 0)
  // 	    {
  // 	      e[recs[i]]=0;
  // 	      p[recs[i]]=0;
  // 	      h[recs[i]]=0;
  // 	    }
  // 	  else
  // 	    {
  // 	      e[recs[i]]=accumulate(eis.begin(),eis.begin()+i,0l);
  // 	      p[recs[i]]=accumulate(pis.begin(),pis.begin()+i,0l);
  // 	      h[recs[i]]=accumulate(his.begin(),his.begin()+i,0l);
  // 	    }
  // 	}
  //   }
  // else
  //   {
  //     for( unsigned i = 0 ; i < recs.size() ; ++i )
  // 	{
  // 	  e[recs[i]]=eis[i];
  // 	  p[recs[i]]=pis[i];
  // 	  h[recs[i]]=his[i];
  // 	}
  //   }
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

size_t simindex::size() const
{
  return e.size();
}
