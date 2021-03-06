#include <Rcpp.h>
#include <zlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <functional>
#include <diseaseSims/ccintermediate.hpp>
#include <SimPopData.hpp>
#include <fwdpp/IO.hpp>

using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]

SimPopData::SimPopData(const std::string & filename,
		       const unsigned long & offset) : mutations(mlist()),
						       gametes(glist()),
						       diploids(dipvector())
{
  this->readPop(filename.c_str(),offset);
}

void SimPopData::readPop(const char * filename,
			 const unsigned long & offset)
{
  gzFile gzin=gzopen(filename,"rb");
  if(gzin==NULL)
    {
      ostringstream error;
      error <<"SimPopData::readPop error: "
	    << filename
	    << " could not be opened for reading\n";
      Rcpp::stop( error.str() );
      return;
    }
  gzseek( gzin,offset,0 );
  KTfwd::read_binary_pop( &(this->gametes), 
			  &(this->mutations),
			  &(this->diploids), std::bind(gzmreader(),std::placeholders::_1),gzin );
  gzclose(gzin);
}

/*
  Returns the positions of neutral mutations on
  the two haplotypes of the ith diploid
*/
List SimPopData::neutralGenotype(const size_t & i) const
{
  if( i == 0 || i > diploids.size() )
    {
      ostringstream error;
      error << "Error: index " << i 
	    << " out of bounds.  Returning empty list\n";
      Rcpp::stop(error.str());
      return List::create();
    }
  std::vector<double> hap1,hap2;
  for(gtype::mcont_const_iterator itr = diploids[i-1].first->mutations.begin() ; 
      itr < diploids[i-1].first->mutations.end() ; ++itr )
    {
      hap1.push_back( (*itr)->pos );
    }

  for(gtype::mcont_const_iterator itr = diploids[i-1].second->mutations.begin() ; 
      itr < diploids[i-1].second->mutations.end() ; ++itr )
    {
      hap2.push_back( (*itr)->pos );
    }
  return List::create( Named("hap1") = hap1,
		       Named("hap2") = hap2 );
}

/*
  Returns the positions of selected mutations on
  the two haplotypes of the ith diploid
*/
List SimPopData::selectedGenotype(const size_t & i) const
{
  if( i == 0 || i > diploids.size() )
    {
      ostringstream error;
      error << "Error: index " << i 
	    << " out of bounds.  Returning empty list\n";
      Rcpp::stop( error.str() );
      return List::create();
    }
  std::vector<double> hap1,hap2;
  for(gtype::mcont_const_iterator itr = diploids[i-1].first->smutations.begin() ; 
      itr < diploids[i-1].first->smutations.end() ; ++itr )
    {
      hap1.push_back( (*itr)->pos );
    }

  for(gtype::mcont_const_iterator itr = diploids[i-1].second->smutations.begin() ; 
      itr < diploids[i-1].second->smutations.end() ; ++itr )
    {
      hap2.push_back( (*itr)->pos );
    }
  return List::create( Named("hap1") = hap1,
		       Named("hap2") = hap2 );
}

dipvector::size_type SimPopData::popsize() const
{
  return diploids.size();
}

RCPP_MODULE(SimPopData)
{
  class_<SimPopData>( "SimPopData" )
    .constructor<string,unsigned long>("Constructor takes a file name and an offset (in bytes")
    .method("ngeno",&SimPopData::neutralGenotype,"Returns a list containing the positions of neutral markers on each haplotype of the i-th diploid.  The argument must be an integer in the range 1 <= x <= population size.")
    .method("sgeno",&SimPopData::selectedGenotype,"Returns a list containing the positions of selected markers on each haplotype of the i-th diploid.  The argument must be an integer in the range 1 <= x <= population size.")
    .method("popsize",&SimPopData::popsize,"Returns the number of diploids.")
    ;
}

//' Read effect sizes from a file at a specific position
//' @param filename The file name.  Should be binary, and either uncompressed or gzip compressed.
//' @param offset The size in bytes where the desired record begins
// [[Rcpp::export]]		 
DataFrame getEsizes( const char * filename,
		     const unsigned long & offset )
{
  NumericVector pos,esize,count,age;
  double row[4];
  unsigned nmuts;

  //PS, zlib is fookin' magic and can read uncompressed files too...
  gzFile gzin = gzopen(filename,"rb");
  if( gzin == NULL )
    {
      ostringstream error;
      error <<"getEsizes error: "
	    << filename 
	    << " could not be opened for reading\n";
      Rcpp::stop( error.str() );
      return DataFrame::create();
    }

  gzseek( gzin, offset, 0 );
  gzread( gzin, &nmuts, sizeof(unsigned) );

  for( unsigned i = 0 ; i < nmuts ; ++i )
    {
      gzread( gzin,&row[0],4*sizeof(double) );
      pos.push_back(row[0]);
      esize.push_back(row[1]);
      count.push_back(row[2]);
      age.push_back(row[3]);
    }
  gzclose(gzin);

  return DataFrame::create( Named("pos") = pos,
			    Named("esize") = esize,
			    Named("count") = count,
			    Named("age") = age );
}

//' Sample a case/control panel from a population.
//' This is used for "on the fly" analysis in place of the make_case_control program,
//' for cases where the user may not want to write a case/control panel to file.
//' @param popfilename  The name of the population file (gzipped, binary)
//' @param offset  The offset of the record in popfilename
//' @param phenofilename The name of the file recording phenotypes (gzipped, binary)
//' @param phenooffset The offset of the record in phenofilename
//' @param ncontrols The number of controls to sample
//' @param ncases The number of cases to sample
//' @param case_proportion The incidence of the disease.  Thus, an individual whose trait value is >= the 1 - case_proportion'th quantile of trait values in the population is a potential case
//' @param control_range A putative control individual is defined as mean +/- control_range*sd, where mean and sd refer to the distribution of the phenotype in the entire population.
//' @param seed Random number seed
// [[Rcpp::export]]
List sampleCCfromPop( const char * popfilename,
		      const unsigned long & offset,
		      const char * phenofilename,
		      const unsigned long & phenooffset,
		      const unsigned & ncontrols,
		      const unsigned & ncases,
		      const double & case_proportion,
		      const double & control_range,
		      const unsigned & seed)
{
  SimPopData spd(popfilename,offset);
  NumericMatrix phenos = getPheno(phenofilename,phenooffset);
  vector<pair<double,double> > phenos2; //converted for use by external function
  for ( int i = 0 ; i < phenos.nrow() ; ++i )
    {
      phenos2.push_back( make_pair( phenos(i,0),phenos(i,1) ) );
    }

  double cutoff;
  pair<double,double> mean_sd = phenosums(phenos2,case_proportion,&cutoff);
  vector<unsigned> put_controls,put_cases;
  grab_putative_CC(mean_sd,phenos2,control_range,cutoff,put_controls,put_cases);
  shuffle( put_controls.begin(), put_controls.end(), default_random_engine(seed) );
  shuffle( put_cases.begin(), put_cases.end(), default_random_engine(seed) );
  cc_intermediate ccblocks(process_population(spd.diploids,phenos2,
					      put_controls,
					      put_cases,
					      ncontrols,
					      ncases) );

  //Fill up the positions
  vector<double> pos;
  copy( ccblocks.neutral.pbegin(),
	ccblocks.neutral.pend(),
	back_inserter(pos) );
  copy( ccblocks.causative.pbegin(),
	ccblocks.causative.pend(),
	back_inserter(pos) );
  IntegerMatrix genos(ncontrols+ncases,ccblocks.neutral.numsites()+ccblocks.causative.size()),
    burdens(ncontrols+ncases,2);
  //Add the neutral sites into genos
  unsigned I = 0;
  for( unsigned i = 0 ; i < ccblocks.neutral.size() ; i+= 2,++I )
    {
      for( unsigned j = 0 ; j < ccblocks.neutral.numsites() ; ++j )
	{
	  genos(I,j) += (ccblocks.neutral[i][j]=='1') ? 1 : 0;
	  genos(I,j) += (ccblocks.neutral[i+1][j]=='1') ? 1 : 0;
	}
    }
  I = 0;
  //And the causative
  for( unsigned i = 0 ; i < ccblocks.causative.size() ; i+= 2,++I )
    {
      for( unsigned j = 0 ; j < ccblocks.causative.numsites() ; ++j )
	{
	  genos(I,j + ccblocks.neutral.numsites()) += (ccblocks.causative[i][j]=='1') ? 1 : 0;
	  genos(I,j + ccblocks.neutral.numsites()) += (ccblocks.causative[i+1][j]=='1') ? 1 : 0;
	  burdens(I,0) = std::count_if( ccblocks.causative[i].begin(),ccblocks.causative[i].end(),[]( const char & c ) { return c == '1'; } );
	  burdens(I,1) = std::count_if( ccblocks.causative[i+1].begin(),ccblocks.causative[i+1].end(),[]( const char & c ) { return c == '1'; } );
	}
    }
  return List::create( Named("pos") = pos,
		       Named("genos") = genos,
		       Named("burdens") = burdens,
		       Named("phenos") = phenos,
		       Named("ncontrols") = ncontrols,
		       Named("ncases") = ncases,
		       Named("neutral") = ccblocks.neutral.numsites(),
		       Named("causative") = ccblocks.causative.numsites());
}

//' Read case/control panel from a file at a specific position
//' @param filename The file name.  Should be binary, and either uncompressed or gzip compressed.
//' @param offset The size in bytes where the desired record begins
// [[Rcpp::export]]	      
List getCCblock( const char * filename,
		 const unsigned long & offset )
{
  gzFile gzin = gzopen(filename,"rb");
  if( gzin == NULL )
    {
      ostringstream error;
      error <<"getCCblock error: "
	    << filename 
	    << " could not be opened for reading\n";
      Rcpp::stop( error.str() );
      return List::create();
    }
  gzseek( gzin,offset,0 );

  unsigned n[4];
  gzread(gzin,&n[0],4*sizeof(unsigned));

  //Read in mutation positions
  NumericVector pos(n[2]+n[3]);
  gzread(gzin,&pos[0],(n[2]+n[3])*sizeof(double));//holy crap this compiles!?!?!?

  //Matrices are allocated rows,cols
  IntegerMatrix genos(n[0]+n[1],n[2]+n[3]),burdens(n[0]+n[1],2);
  NumericMatrix phenos(n[0]+n[1],2);

  //read in the genotype matrix
  unsigned nones,ntwos;
  for( unsigned i = 0 ; i < n[0]+n[1] ; ++i )
    {
      gzread( gzin,&nones,sizeof(unsigned) );
      vector<unsigned> ones(nones);
      gzread( gzin,&ones[0],nones*sizeof(unsigned) );

      for( unsigned j = 0 ; j < nones ; ++j )
	{
	  genos(i,ones[j]) = 1;
	}

      gzread( gzin,&ntwos,sizeof(unsigned) );
      vector<unsigned> twos(ntwos);
      gzread( gzin,&twos[0],ntwos*sizeof(unsigned) );

      for( unsigned j = 0 ; j < ntwos ; ++j )
	{
	  genos(i,twos[j]) = 2;
	}
    }

  //read in burdens
  unsigned b[2];
  for( unsigned i = 0 ; i < n[0]+n[1] ; ++i )
    {
      gzread(gzin,&b[0],2*sizeof(unsigned));
      burdens(i,0)=b[0];
      burdens(i,1)=b[1];
    }

  //And phenos
  double p[2];
  for( unsigned i = 0 ; i < n[0]+n[1] ; ++i )
    {
      gzread(gzin,&p[0],2*sizeof(double));
      phenos(i,0)=p[0];
      phenos(i,1)=p[1];
    }

  gzclose(gzin);

  //Apparently, Rcpp assigns to columns first...
  //gzread( gzin,&burdens[0],2*(n[0]+n[1])*sizeof(unsigned) );
  //gzread( gzin,&phenos[0],2*(n[0]+n[1])*sizeof(double) );

  return List::create( Named("pos") = pos,
		       Named("genos") = genos,
		       Named("burdens") = burdens,
		       Named("phenos") = phenos,
		       Named("ncontrols") = n[0],
		       Named("ncases") = n[1],
		       Named("neutral") = n[2],
		       Named("causative") = n[3]);
}

//' Read case/control ids from a file at a specific position
//' @param filename The file name.  Should be binary, and either uncompressed or gzip compressed.
//' @param offset The size in bytes where the desired record begins
//' @return A list of who the controls and cases were.  These are values between 1 and N, the population size.
// [[Rcpp::export]]
List getCCids( const char * filename,
	       const unsigned long & offset )
{
  gzFile gzin = gzopen(filename,"rb");
  if( gzin == NULL )
    {
      ostringstream error;
      error <<"getCCids error: "
	    << filename 
	    << " could not be opened for reading\n";
      Rcpp::stop( error.str() );
      return List::create();
    }
  gzseek( gzin,offset,0 );

  unsigned ncontrols,ncases;
  gzread(gzin,&ncontrols,sizeof(unsigned));
  gzread(gzin,&ncases,sizeof(unsigned));

  NumericVector controls(ncontrols),cases(ncases);
  gzread( gzin, &controls[0], ncontrols*sizeof(unsigned) );
  gzread( gzin, &cases[0], ncases*sizeof(unsigned) );

  gzclose( gzin );

  return( List::create( Named("controls") = controls,
			Named("case") = cases ) );
}

//' Read case/control panel from a file at a specific position
//' @param filename The file name.  Should be binary, and either uncompressed or gzip compressed.
//' @param offset The size in bytes where the desired record begins
// [[Rcpp::export]]
NumericMatrix getPheno(const char * filename,
		       const unsigned long & offset)
{
  gzFile gzin = gzopen(filename,"rb");
  if( gzin == NULL )
    {
      ostringstream error;
      error <<"getPheno error: "
	    << filename 
	    << " could not be opened for reading\n";
      Rcpp::stop( error.str() );
      return NumericMatrix();
    }
  gzseek( gzin,offset,0 );
  unsigned n;
  gzread(gzin,&n,sizeof(unsigned));
  
  NumericMatrix rv(n,2);
  double p[2];
  for(unsigned i=0;i<n;++i)
    {
      gzread(gzin,&p[0],2*sizeof(double));
      rv(i,0)=p[0];
      rv(i,1)=p[1];
    }
  gzclose(gzin);
  return rv;
}

