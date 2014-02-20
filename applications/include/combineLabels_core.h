#ifndef combineLabels_core_h
#define combineLabels_core_h

#include <map>
#include <stdexcept>

extern long int ran2Seed;
extern long int ran2initialSeed;

typedef std::map <short, short> countMap;

short int getMostPopular (const countMap& cMap);
bool 	  isEquivocal 	 (const countMap& cMap);
short int decideOnTie 	 (const countMap& cMap) throw (std::invalid_argument);

#endif // combineLabels_core_h
