#include <diseaseSims/util.hpp>
#include <fwdpp/IO.hpp>

using namespace std;

popstruct readPop( gzFile gzin )
{
  popstruct p;
  KTfwd::read_binary_pop( &p.gametes, &p.mutations, &p.diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
  return p;
}

vector<double> getG( const dipvector & diploids,
		     const Gfxn_t & dipG )
{
  vector<double> rv;
  for_each( diploids.begin(),diploids.end(),[&rv,&dipG](const diploid_t & __d ) { rv.push_back( dipG(__d.first,__d.second) ); } );
  return rv;
}

vmcount_t get_mut_counts( const glist::const_iterator & g1,
			  const glist::const_iterator & g2,
const bool & selected )
{
  vmcount_t rv;

  auto updater = [&rv](const mlist::iterator & __mut) {
    auto __itr =  find_if(rv.begin(),rv.end(),[&__mut](const pair<mlist::iterator,unsigned> & __p) {
	return __p.first == __mut;
      } );
    if(__itr == rv.end())
      rv.push_back(make_pair(__mut,1u));
    else
      __itr->second++;
  };
  if(selected)
    {
      for_each( g1->smutations.begin(), g1->smutations.end(),updater );
      for_each( g2->smutations.begin(), g2->smutations.end(),updater );
    }
  else
    {
      for_each( g1->mutations.begin(), g1->mutations.end(),updater );
      for_each( g2->mutations.begin(), g2->mutations.end(),updater );
    }
  return rv;
}
