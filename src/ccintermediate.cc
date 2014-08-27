#include <ccintermediate.hpp>
#include <mutation_with_age.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <boost/bind.hpp>

using namespace std;

cc_intermediate::cc_intermediate(void) : neutral(Sequence::SimData()),
					 causative(Sequence::SimData()),
					 min_n( vector<char>() ),
					 min_c( vector<char>() ),
					 phenotypes( vector<pair<double,double> >() )
{
}

//called by process_subset to update data blocks
void update_block( glist::const_iterator::value_type::mutation_container::const_iterator beg,  //gamete1
		   glist::const_iterator::value_type::mutation_container::const_iterator end,
		   glist::const_iterator::value_type::mutation_container::const_iterator beg2, //gamete2
		   glist::const_iterator::value_type::mutation_container::const_iterator end2,
		   const unsigned & i,
		   vector< pair<double,string> > & datablock,
		   const unsigned & ttl,
		   const unsigned & offset)
{
  std::function<bool(const std::pair<double,std::string> &, const double &)> sitefinder = [](const std::pair<double,std::string> & site,
											     const double & d ) 
    {
      return std::fabs(site.first-d) <= std::numeric_limits<double>::epsilon();
    };
  for( ; beg < end ; ++beg )
    {
      double mutpos = (*beg)->pos;
      vector< pair<double,string> >::iterator itr =  find_if(datablock.begin(),
							     datablock.end(),
							     std::bind(sitefinder,std::placeholders::_1,mutpos));
      if( itr == datablock.end() )
	{
	  datablock.push_back( make_pair(mutpos,string(ttl,'0')) );
	  datablock[datablock.size()-1].second[offset + 2*i] = '1';
	}
      else
	{
	  assert(offset+2*i < itr->second.size());
	  itr->second[offset + 2*i] = '1';
	}
    }
  for( ; beg2 < end2 ; ++beg2 )
    {
      double mutpos = (*beg2)->pos;
      vector< pair<double,string> >::iterator itr =  find_if(datablock.begin(),
							     datablock.end(),
							     std::bind(sitefinder,std::placeholders::_1,mutpos));
      if( itr == datablock.end() )
	{
	  datablock.push_back( make_pair(mutpos,string(ttl,'0')) );
	  datablock[datablock.size()-1].second[offset + 2*i + 1] = '1';
	}
      else
	{
	  assert(offset+2*i + 1 < itr->second.size());
	  itr->second[offset + 2*i + 1] = '1';
	}
    }
}

void process_subset( vector< pair<double,string> > & datablock_neut,
		     vector< pair<double,string> > & datablock_sel,
		     vector< pair<double,double> > & ccphenos,
		     const vector< pair<glist::iterator,glist::iterator> > & diploids,
		     const vector<pair<double,double> > & popphenos,
		     const vector<size_t> & indlist,
		     const unsigned & maxnum,
		     const unsigned & ttl,
		     const unsigned & offset)
{
  vector< pair<double,string> >::iterator itr;

  for( unsigned i = 0 ; i < maxnum ; ++i )
    {
      ccphenos.push_back( popphenos[ indlist[i] ] );
      update_block(  diploids[ indlist[i] ].first->mutations.begin(),
		     diploids[ indlist[i] ].first->mutations.end(),
		     diploids[ indlist[i] ].second->mutations.begin(),
		     diploids[ indlist[i] ].second->mutations.end(),
		     i,
		     datablock_neut,
		     ttl,offset
		     );
      update_block(  diploids[ indlist[i] ].first->smutations.begin(),
		     diploids[ indlist[i] ].first->smutations.end(),
		     diploids[ indlist[i] ].second->smutations.begin(),
		     diploids[ indlist[i] ].second->smutations.end(),
		     i,
		     datablock_sel,
		     ttl,offset
		     );

      /*
      //Old code replaced by update_block
      //neutral
      for( unsigned mut = 0 ; mut < diploids[ indlist[i] ].first->mutations.size() ; ++mut )
	{
	  double mutpos =  diploids[ indlist[i] ].first->mutations[mut]->pos;
	  itr = find_if(datablock_neut.begin(),
			datablock_neut.end(),
			std::bind(KTfwd::find_mut_pos(),std::placeholders::_1,mutpos));
	  if( itr == datablock_neut.end() )
	    {
	      datablock_neut.push_back( make_pair(mutpos,string(ttl,'0')) );
	      datablock_neut[datablock_neut.size()-1].second[offset + 2*i] = '1';
	    }
	  else
	    {
	      assert(offset+2*i < itr->second.size());
	      itr->second[offset + 2*i] = '1';
	    }
	}
      for( unsigned mut = 0 ; mut < diploids[ indlist[i] ].second->mutations.size() ; ++mut )
	{
	  double mutpos =  diploids[ indlist[i] ].second->mutations[mut]->pos;
	  itr = find_if(datablock_neut.begin(),
			datablock_neut.end(),
			std::bind(KTfwd::find_mut_pos(),std::placeholders::_1,mutpos));
	  if( itr == datablock_neut.end() )
	    {
	      datablock_neut.push_back( make_pair(mutpos,string(ttl,'0')) );
	      datablock_neut[datablock_neut.size()-1].second[2*i + 1] = '1';
	    }
	  else
	    {
	      assert( (offset + 2*i + 1) < itr->second.size() );
	      itr->second[offset + 2*i + 1] = '1';
	    }
	}
      //selected
      for( unsigned mut = 0 ; mut < diploids[ indlist[i] ].first->smutations.size() ; ++mut )
	{
	  double mutpos =  diploids[ indlist[i] ].first->smutations[mut]->pos;
	  itr = find_if(datablock_sel.begin(),
			datablock_sel.end(),
			std::bind(KTfwd::find_mut_pos(),std::placeholders::_1,mutpos));
	  if( itr == datablock_sel.end() )
	    {
	      datablock_sel.push_back( make_pair(mutpos,string(ttl,'0')) );
	      datablock_sel[datablock_sel.size()-1].second[offset + 2*i] = '1';
	    }
	  else
	    {
	      assert( offset+2*i < itr->second.size() );
	      itr->second[offset+2*i] = '1';
	    }
	}
      for( unsigned mut = 0 ; mut < diploids[ indlist[i] ].second->smutations.size() ; ++mut )
	{
	  double mutpos =  diploids[ indlist[i] ].second->smutations[mut]->pos;
	  itr = find_if(datablock_sel.begin(),
			datablock_sel.end(),
			std::bind(KTfwd::find_mut_pos(),std::placeholders::_1,mutpos));
	  if( itr == datablock_sel.end() )
	    {
	      datablock_sel.push_back( make_pair(mutpos,string(ttl,'0')) );
	      datablock_sel[datablock_sel.size()-1].second[2*i + 1] = '1';
	    }
	  else
	    {
	      assert( (offset + 2*i + 1) < itr->second.size() );
	      itr->second[offset + 2*i + 1] = '1';
	    }
	}
      */
    }
}

cc_intermediate process_population( const vector< pair<glist::iterator,glist::iterator> > & diploids,
				    const vector<pair<double,double> > & phenotypes,
				    const vector<size_t> & put_controls,
				    const vector<size_t> & put_cases,
				    const unsigned & ncontrols,
				    const unsigned & ncases)
{
  cc_intermediate rv;

  /*
    vector of position,genotype pairs. 
    For each position, there will be ncontrols + 
    ncases genotypes.  Controls before cases.

    neutral = neutral mutations in population
    selected = causative mutation in population
    We keep the 2 classes separate for easier
    processing downstream.
  */
  vector< pair<double,string> > neutral,selected;

  //Go thru controls first
  process_subset( neutral, selected,
		  rv.phenotypes,
		  diploids,
		  phenotypes,
		  put_controls,
		  ncontrols,
		  2*(ncontrols+ncases),
		  0 );
  //cases
  process_subset( neutral, selected,
		  rv.phenotypes,
		  diploids,
		  phenotypes,
		  put_cases,
		  ncases,
		  2*(ncontrols+ncases),
		  2*ncontrols);

  sort( neutral.begin(), neutral.end(), 
		[](std::pair<double,std::string> lhs,
		   std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
  //std::bind(KTfwd::sortpos(),std::placeholders::_1,std::placeholders::_2) );
  sort( selected.begin(), selected.end(), 
	[](std::pair<double,std::string> lhs,
	   std::pair<double,std::string> rhs) { return lhs.first < rhs.first; });
  //std::bind(KTfwd::sortpos(),std::placeholders::_1,std::placeholders::_2) );

  rv.neutral.assign( neutral.begin(), neutral.end() );  
  rv.causative.assign( selected.begin(), selected.end() );  

  Sequence::RemoveInvariantColumns(&rv.neutral);
  Sequence::RemoveInvariantColumns(&rv.causative);

  //Define the minor allele state

  //The block below deterimines minor allele based on controls
  /*
  for( Sequence::SimData::const_site_iterator i = rv.neutral.sbegin() ; 
       i < rv.neutral.send() ; ++i )
    {
      size_t c = count(i->second.begin(),i->second.begin() + 2*ncontrols,'1');
      rv.min_n.push_back( (c < ncontrols) ? '1' : '0' );
    }
  for( Sequence::SimData::const_site_iterator i = rv.causative.sbegin() ; 
       i < rv.causative.send() ; ++i )
    {
      size_t c = count(i->second.begin(),i->second.begin() + 2*ncontrols,'1');
      rv.min_c.push_back( (c < ncontrols) ? '1' : '0' );
    }
  */

  /*
    This block is based on whole case/control sample.
    This is actually what was done in Thornton, Foran, and Long (2013)
  */
  for( Sequence::SimData::const_site_iterator i = rv.neutral.sbegin() ; 
       i < rv.neutral.send() ; ++i )
    {
      size_t c = count(i->second.begin(),i->second.end(),'1');
      rv.min_n.push_back( (c <= rv.neutral.size()/2) ? '1' : '0' );
    }
  for( Sequence::SimData::const_site_iterator i = rv.causative.sbegin() ; 
       i < rv.causative.send() ; ++i )
    {
      size_t c = count(i->second.begin(),i->second.end(),'1');
      rv.min_c.push_back( (c < rv.causative.size()/2) ? '1' : '0' );
    }
  return rv;
}
