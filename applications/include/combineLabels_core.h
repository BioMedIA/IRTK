#ifndef combineLabels_core_h
#define combineLabels_core_h

#include <map>
#include <stdexcept>

extern char * output_name;
extern char * * input_names;
extern char * mask_name;
extern long int ran2Seed;
extern long int ran2initialSeed;

typedef std::map <short, short> countMap;
extern 	std::map <short, short>::iterator iter;

short int getMostPopular (countMap cMap);
bool 	  isEquivocal 	 (countMap cMap);
short int decideOnTie 	 (countMap cMap) throw (std::invalid_argument);

#endif // combineLabels_core_h
