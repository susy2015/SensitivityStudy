#ifndef _SSTrimCommon_H_
#define _SSTrimCommon_H_

#include <string>
#include <iostream>
#include <vector>

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/baselineDef.h"

inline size_t find_Nth
(
  const std::string & str ,   // where to work
  unsigned            N ,     // N'th ocurrence
  const std::string & find    // what to 'find'
)
{
  if ( 0==N ) { return std::string::npos; }
  size_t pos,from=0;
  unsigned i=0;
  while ( i<N )
  {
    pos=str.find(find,from);
    if ( std::string::npos == pos ) { break; }
    from = pos + 1; // from = pos + find.size();
    ++i;
  }
  return pos;
}

inline double GetHTTops( std::vector<TLorentzVector> vTops )
{
  int HTTops = 0;

  for(auto Top : vTops)
  {
    HTTops += Top.Pt();
  }
  return HTTops;
}       // -----  end of function VarPerEvent::GetnTops  -----


#endif
