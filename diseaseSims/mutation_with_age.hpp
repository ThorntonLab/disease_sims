#ifndef __MUTATION_WITH_AGE_HPP__
#define __MUTATION_WITH_AGE_HPP__

#include <fwdpp/forward_types.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <cmath>
#include <iosfwd>
#include <sstream>

#include <zlib.h>

#ifndef USE_STANDARD_CONTAINERS
#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>
#include <boost/pool/pool_alloc.hpp>
#include <boost/unordered_set.hpp>
#include <boost/functional/hash.hpp>
#else
#include <vector>
#include <list>
#include <unordered_set>
#include <functional>
#endif

struct mutation_with_age : public KTfwd::mutation_base
{
  mutable unsigned o;
  double s;
  char label; //'A' = "amino acid", 'S' = "synonymous"  Only used here in that A = causative, S = neutral.
  mutation_with_age( const double & position, const double & sel_coeff,
		     const unsigned & count,const unsigned & origin, const char & ch,
		     const bool & n=true) 
    : mutation_base(position,count,n),o(origin),s(sel_coeff),label(ch)
  {
  }
  bool operator==(const mutation_with_age & rhs) const
  {
    return( std::fabs(this->pos-rhs.pos) <= std::numeric_limits<double>::epsilon() &&
	    this->s == rhs.s );
  }	
};

typedef mutation_with_age TFLmtype;
//boost containers
#ifndef USE_STANDARD_CONTAINERS
typedef boost::pool_allocator<TFLmtype> mut_allocator;
typedef boost::container::list<TFLmtype,mut_allocator > mlist;
typedef KTfwd::gamete_base<TFLmtype,mlist> gtype;
typedef boost::pool_allocator<gtype> gam_allocator;
typedef boost::container::list<gtype,gam_allocator > glist;
typedef boost::container::vector<TFLmtype> mvector;
typedef boost::container::vector<unsigned> ftvector;
typedef boost::unordered_set<double,boost::hash<double>,KTfwd::equal_eps > lookup_table_type;
#else
typedef std::list<TFLmtype > mlist;
typedef KTfwd::gamete_base<TFLmtype, mlist> gtype;
typedef std::list<gtype> glist;
typedef std::vector<TFLmtype> mvector;
typedef std::vector<unsigned> ftvector;
typedef std::unordered_set<double,std::hash<double>,KTfwd::equal_eps > lookup_table_type;
#endif

//function object to write mutation data in binary format
struct mwriter
{
  typedef void result_type;
  result_type operator()( const TFLmtype & m, std::ostringstream & buffer ) const
  {
    buffer.write( reinterpret_cast< const char * >(&m.n),sizeof(unsigned) );
    buffer.write( reinterpret_cast< const char * >(&m.o),sizeof(unsigned) );
    buffer.write( reinterpret_cast< const char * >(&m.neutral),sizeof(bool) );
    buffer.write( reinterpret_cast< const char * >(&m.pos),sizeof(double) );
    buffer.write( reinterpret_cast< const char * >(&m.s),sizeof(double) );
    buffer.write( reinterpret_cast< const char * >(&m.label),sizeof(char) );
  }
};

//function object to read mutation data in binary format
struct mreader
{
  typedef TFLmtype result_type;
  result_type operator()( std::istream & in ) const
  {
    unsigned n;
    in.read( reinterpret_cast< char * >(&n),sizeof(unsigned) );
    unsigned o;
    in.read( reinterpret_cast< char * >(&o),sizeof(unsigned) );
    bool neut;
    in.read( reinterpret_cast< char * >(&neut),sizeof(bool) );
    double pos;
    in.read( reinterpret_cast< char * >(&pos),sizeof(double) );
    double s;
    in.read( reinterpret_cast< char * >(&s),sizeof(double) );
    char label;
    in.read( reinterpret_cast< char * >(&label),sizeof(char) );
    return result_type(pos,s,n,o,label,neut);
  }
};

//function object to read mutation data in binary format from a gzipped file
struct gzmreader
{
  typedef TFLmtype result_type;
  result_type operator()( gzFile gzin ) const
  {
    unsigned n;
    gzread(gzin,&n,sizeof(unsigned));
    unsigned o;
    gzread(gzin,&o,sizeof(unsigned));
    bool neut;
    gzread(gzin,&neut,sizeof(bool));
    double pos;
    gzread(gzin,&pos,sizeof(double));
    double s;
    gzread(gzin,&s,sizeof(double));
    char label;
    gzread(gzin,&label,sizeof(char));
    return result_type(pos,s,n,o,label,neut);
  }
};

#endif
