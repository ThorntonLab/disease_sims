#include <diseaseSims/ccintermediate.hpp>
#include <diseaseSims/mutation_with_age.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <fwdpp/sampling_functions.hpp>
#include <cassert>

#include <gsl/gsl_statistics.h>

using namespace std;

cc_intermediate::cc_intermediate(void) : ncontrols(0),ncases(0),
					 neutral(Sequence::SimData()),
					 causative(Sequence::SimData()),
					 min_n( vector<char>() ),
					 min_c( vector<char>() ),
					 G ( vector<double>() ),
					 E ( vector<double>() ),
					 control_ids( vector<unsigned>() ),
					 case_ids( vector<unsigned>() )
{
}

std::ostream & cc_intermediate::buffer( std::ostream & ccbuffer ) const
{
  ccbuffer.write( reinterpret_cast<const char *>(&ncontrols),sizeof(unsigned) );
  ccbuffer.write( reinterpret_cast<const char *>(&ncases),sizeof(unsigned) );
  unsigned temp = neutral.numsites();
  ccbuffer.write( reinterpret_cast<char *>(&temp),sizeof(unsigned) );
  temp = causative.numsites();
  ccbuffer.write( reinterpret_cast<char *>(&temp),sizeof(unsigned) );

  for( Sequence::SimData::const_pos_iterator p = neutral.pbegin() ; 
       p < neutral.pend() ; ++p )
    {
      double x = *p;
      ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
    }

  for( Sequence::SimData::const_pos_iterator p = causative.pbegin() ; 
       p < causative.pend() ; ++p )
    {
      double x = *p;
      ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
    }

  //iterate over the diploids and write block for association tests
  for( unsigned ind = 0 ; ind < neutral.size() ; ind += 2 ) 
    {
      vector<unsigned> ones,twos;
      //neutral genotypes for this individual
      for( unsigned site = 0 ; site < neutral.numsites() ; ++site )
	{
	  //count # copies of minor allele at this site in this individual
	  unsigned cminor = ( (neutral[ind][site] == min_n[site]) ? 1 : 0 ) +
	    ( (neutral[ind+1][site] == min_n[site]) ? 1 : 0 );
	  if( cminor==1 )
	    {
	      ones.push_back(site);
	    }
	  else if ( cminor == 2 )
	    {
	      twos.push_back(site);
	    }
	  //ccbuffer.write( reinterpret_cast<char *>(&cminor), sizeof(unsigned) );
	}
      //causative genotypes for this individual
      for( unsigned site = 0 ; site < causative.numsites() ; ++site )
	{
	  //count # copies of minor allele at this site in this individual
	  unsigned cminor = ( (causative[ind][site] == min_c[site]) ? 1 : 0 ) +
	    ( (causative[ind+1][site] == min_c[site]) ? 1 : 0 );
	  //ccbuffer.write( reinterpret_cast<char *>(&cminor), sizeof(unsigned) );
	  if( cminor==1 )
	    {
	      ones.push_back(neutral.numsites()+site);
	    }
	  else if ( cminor == 2 )
	    {
	      twos.push_back(neutral.numsites()+site);
	    }
	}
      //update the buffer
      unsigned n = ones.size();
      ccbuffer.write( reinterpret_cast<char *>(&n),sizeof(unsigned) );
      ccbuffer.write( reinterpret_cast<char *>(&ones[0]),n*sizeof(unsigned) );
      n = twos.size();
      ccbuffer.write( reinterpret_cast<char *>(&n),sizeof(unsigned) );
      ccbuffer.write( reinterpret_cast<char *>(&twos[0]),n*sizeof(unsigned) );
    }

  //Now, output # of causative mutations on each haplotype carried by this diploid
  for( unsigned ind = 0 ; ind < causative.size() ; ind += 2 ) 
    {
      unsigned ncaus = count( causative[ind].begin(),
			      causative[ind].end(),'1' );
      ccbuffer.write( reinterpret_cast<char *>(&ncaus), sizeof(unsigned) );
      ncaus = count( causative[ind+1].begin(),
		     causative[ind+1].end(),'1' );
      ccbuffer.write( reinterpret_cast<char *>(&ncaus), sizeof(unsigned) );
    }
  //Give phenos of controls & cases
  for(unsigned i = 0 ; i < G.size() ; ++i )
    {
      ccbuffer.write(reinterpret_cast<const char *>(&G[i]),sizeof(double));
      ccbuffer.write(reinterpret_cast<const char *>(&E[i]),sizeof(double));
    }
  // for( vector< pair<double,double> >::const_iterator i = phenotypes.begin() ; 
  //      i != phenotypes.end() ; ++i )
  //   {
  //     double x = i->first;
  //     ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
  //     x = i->second;
  //     ccbuffer.write( reinterpret_cast<char *>(&x),sizeof(double) );
  //   }
  return ccbuffer;
}

std::ostream & operator<<(std::ostream & o, const cc_intermediate & b )
{
  return b.buffer(o);
}

//called by process_subset to update data blocks
void update_block( poptype::gamete_t::mutation_container::const_iterator beg,  //gamete1
		   poptype::gamete_t::mutation_container::const_iterator end,
		   poptype::gamete_t::mutation_container::const_iterator beg2, //gamete2
		   poptype::gamete_t::mutation_container::const_iterator end2,
		   const poptype & pop,
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
      double mutpos = pop.mutations[*beg].pos;
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
      double mutpos = pop.mutations[*beg2].pos;
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
		     //vector< pair<double,double> > &cphenos,
		     vector<double> & ccG,
		     vector<double> & ccE,
		     const poptype & pop,
		     const vector<pair<double,double> > & popphenos,
		     const vector<unsigned> & indlist,
		     const unsigned & maxnum,
		     const unsigned & ttl,
		     const unsigned & offset,
		     vector<unsigned> & IDS)
{
  vector< pair<double,string> >::iterator itr;
  for( unsigned i = 0 ; i < maxnum ; ++i )
    {
      assert( i < maxnum );
      IDS.push_back(indlist[i]); //record the individual who is a case or control
      assert(indlist[i] < diploids.size());
      //ccphenos.push_back( popphenos[ indlist[i] ] );
      ccG.push_back( popphenos[ indlist[i] ].first );
      ccE.push_back( popphenos[ indlist[i] ].second );
      update_block(  pop.gametes[pop.diploids[ indlist[i] ].first].mutations.begin(),
		     pop.gametes[pop.diploids[ indlist[i] ].first].mutations.end(),
		     pop.gametes[pop.diploids[ indlist[i] ].second].mutations.begin(),
		     pop.gametes[pop.diploids[ indlist[i] ].second].mutations.end(),
		     pop,
		     i,
		     datablock_neut,
		     ttl,offset
		     );
      update_block(  pop.gametes[pop.diploids[ indlist[i] ].first].smutations.begin(),
		     pop.gametes[pop.diploids[ indlist[i] ].first].smutations.end(),
		     pop.gametes[pop.diploids[ indlist[i] ].second].smutations.begin(),
		     pop.gametes[pop.diploids[ indlist[i] ].second].smutations.end(),
		     pop,
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

cc_intermediate process_population( const poptype & pop,
				    const vector<pair<double,double> > & phenotypes,
				    const vector<unsigned> & put_controls,
				    const vector<unsigned> & put_cases,
				    const unsigned & ncontrols,
				    const unsigned & ncases)
{
  cc_intermediate rv;
  rv.ncontrols = ncontrols;
  rv.ncases = ncases;

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
		  rv.G,rv.E,
		  pop,
		  phenotypes,
		  put_controls,
		  ncontrols,
		  2*(ncontrols+ncases),
		  0, rv.control_ids );
  //cases
  process_subset( neutral, selected,
		  rv.G,rv.E,
		  pop,
		  phenotypes,
		  put_cases,
		  ncases,
		  2*(ncontrols+ncases),
		  2*ncontrols, rv.case_ids);

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

  Sequence::removeInvariantPos(rv.neutral);
  Sequence::removeInvariantPos(rv.causative);

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

std::pair<double,double> phenosums(const vector<pair<double,double> > & phenos, const double & case_proportion, double * cutoff)
{
  vector<double> pcopy;
  for(unsigned i=0;i<phenos.size();++i)
    {
      pcopy.push_back(phenos[i].first+phenos[i].second);
    }
  sort(pcopy.begin(),pcopy.end());
  //get the upper case_proportion-th'd quantile from the pheno dist
  *cutoff = gsl_stats_quantile_from_sorted_data(&pcopy[0],1,pcopy.size(),1.-case_proportion);
  return std::make_pair( gsl_stats_mean(&pcopy[0],1,pcopy.size()),
			 gsl_stats_sd(&pcopy[0],1,pcopy.size()) );
}

void grab_putative_CC( const pair<double,double> & mean_sd,
		       const vector<pair<double,double> > & phenotypes,
		       const double & crange,
		       const double & cutoff,
		       std::vector<unsigned> & put_controls,
		       std::vector<unsigned> & put_cases )
{
  put_controls.clear();
  put_cases.clear();
  for( unsigned i = 0 ; i < phenotypes.size() ; ++i )
    {
      const double P = phenotypes[i].first+phenotypes[i].second;
      if( P >= cutoff )
	{
	  put_cases.push_back(i);
	}
      /*
	Issue alert!!!
	TFL (2013) claim that "controls are within 1 standard
	deviation of the mean", by which they mean in terms of phenotype.
	
	Well, that was technically true, but a problem for reproducibility
	because our code to make case/control panels required that a putative
	control's phenotype be within 0.5*sd of population mean phenotype!!!
      */
      else if ( P >= mean_sd.first - crange*mean_sd.second &&
		P <= mean_sd.first + crange*mean_sd.second )
	{
	  put_controls.push_back(i);
	}
    }

  sort( put_controls.begin(), put_controls.end() );
  sort( put_cases.begin(), put_cases.end() );
}

std::vector<int> decode(const std::vector<std::int8_t> & d)
{
  std::vector<int> rv;
  for( unsigned pos = 0 ; pos < d.size() ; ++pos )
    {
      for(unsigned byte=0;byte<8;byte+=2)
        {
          bool i = (d[pos]&(1<<byte));
          bool j = (d[pos]&(1<<(byte+1)));
          if(!i && !j) return rv;
          else if (!i && j ) rv.push_back(0);
          else if (i && !j ) rv.push_back(1);
          else if (i&&j) rv.push_back(2);
        }
    }
  return rv;
}

std::vector<std::int8_t> encode( const std::vector<std::string> & neutral,
				 const std::vector<std::string> & causative)
{
  std::vector<std::int8_t> rv;
  unsigned nsites = (neutral.empty()) ? 0 : neutral[0].size(),
    csites = (causative.empty()) ? 0 : causative[0].size();
  
  for( unsigned ind = 0 ; ind < neutral.size() ; ind += 2 )
    {
      std::int8_t i8=0;
      unsigned byte = 0;
      for( unsigned site = 0 ; site < nsites ;++site )
        {
          unsigned c = (neutral[ind][site] == '1' ? 1 : 0) +  (neutral[ind+1][site] == '1' ? 1 : 0);
          if( c== 0)
            {
              i8 |= (1 <<(byte+1));
            }
          else if ( c == 1 )
            {
              i8 |= (1 <<(byte));
            }
          else if (c==2)
            {
              i8 |= (1 <<(byte));
              i8 |= (1 <<(byte+1));
            }
          byte+=2;
          if(byte==8 || (site == nsites-1 && causative.empty()))
            {
              rv.push_back(i8);
              i8=0;
              byte=0;
            }
        }
      for( unsigned site = 0 ; site < csites ; ++site )
	{
	  unsigned c = (causative[ind][site] == '1' ? 1 : 0) +  (causative[ind+1][site] == '1' ? 1 : 0);
          if( c== 0)
            {
              i8 |= (1 <<(byte+1));

            }
          else if ( c == 1 )
            {
              i8 |= (1 <<(byte));
            }
          else if (c==2)
            {
              i8 |= (1 <<(byte));
              i8 |= (1 <<(byte+1));
            }
          byte+=2;
          if(byte==8 || site == csites-1)
            {
              rv.push_back(i8);
              i8=0;
              byte=0;
            }
	}
    }
  return rv;
}

