#ifndef __SIMINDEX_HPP__
#define __SIMINDEX_HPP__

#include <map>
#include <vector>
#include <cstddef>

class simindex
{
private:
  std::map<unsigned,long> e,p,h;
  bool fileproblem;

  bool mono_increasing( const std::vector<long> & v ) const;
public:
  simindex(const char * filename);

  bool file_problem() const;
  bool eexists(const unsigned & i) const;
  bool pexists(const unsigned & i) const;
  bool hexists(const unsigned & i) const;
  long eoffset(const unsigned & i) const;
  long poffset(const unsigned & i) const;
  long hoffset(const unsigned & i) const;

  //std::map<unsigned,long>::size_type size() const;
  size_t size() const;
};

#endif
