#include <diseaseSims/util.hpp>
#include <fwdpp/sugar/serialization.hpp>

using namespace std;

poptype readPop( gzFile gzin )
{
  poptype p(0);
  //KTfwd::read_binary_pop( &p.gametes, &p.mutations, &p.diploids, std::bind(gzmreader(),std::placeholders::_1),gzin );
  KTfwd::gzdeserialize()(gzin,p,std::bind(KTfwd::mutation_reader<TFLmtype>(),std::placeholders::_1));
  return p;
}

vector<double> getG( const poptype & pop,
		     const Gfxn_t & dipG )
{
  vector<double> rv;
  for(const auto & dip : pop.diploids)
    {
      rv.push_back(dipG(pop.gametes[dip.first],pop.gametes[dip.second],pop.mutations));
    }
  return rv;
}

vmcount_t get_mut_counts( const poptype & pop,
			  const std::size_t dipindex,
			  const bool & selected )
{
  vmcount_t rv;

  auto updater = [&rv](const std::size_t __mut) {
    auto __itr =  find_if(rv.begin(),rv.end(),[&__mut](const pair<std::size_t,unsigned> & __p) {
	return __p.first == __mut;
      } );
    if(__itr == rv.end())
      rv.push_back(make_pair(__mut,1u));
    else
      __itr->second++;
  };
  if(selected)
    {
      for_each( pop.gametes[pop.diploids[dipindex].first].smutations.begin(),
		pop.gametes[pop.diploids[dipindex].first].smutations.end(),
		updater );
      for_each( pop.gametes[pop.diploids[dipindex].second].smutations.begin(),
		pop.gametes[pop.diploids[dipindex].second].smutations.end(),
		updater );
    }
  else
    {
      for_each( pop.gametes[pop.diploids[dipindex].first].mutations.begin(),
		pop.gametes[pop.diploids[dipindex].first].mutations.end(),
		updater );
      for_each( pop.gametes[pop.diploids[dipindex].second].mutations.begin(),
		pop.gametes[pop.diploids[dipindex].second].mutations.end(),
		updater );
    }
  return rv;
}
