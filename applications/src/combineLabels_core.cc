// avoid tweaking build system for the moment
#include "../include/combineLabels_core.h"
#include <irtkImage.h>
#include <nr.h>

#include <sys/types.h>
#include <stdexcept>

#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif


long ran2Seed;
long ran2initialSeed;


short getMostPopular(const countMap& cMap){
  short maxCount = 0, mostPopLabel = -1;

  if ( cMap.empty() ) {
    // No votes to count, treat as background.
    return 0;
  }

  countMap::const_iterator iter;
  const countMap::const_iterator endMap = cMap.end();
  for (iter = cMap.begin(); iter != endMap; ++iter){

    if (iter->second > maxCount){
      maxCount     = iter->second;
      mostPopLabel = iter->first;
    }
  }

  return mostPopLabel;
}

bool isEquivocal(const countMap& cMap){
  short maxCount = 0;
  short numberWithMax = 0;

  if ( cMap.empty() ) {
    // No votes to count, treat as background.
    return false;
  }

  countMap::const_iterator iter;
  const countMap::const_iterator endMap = cMap.end();
  for (iter = cMap.begin(); iter != endMap; ++iter){
    if (iter->second > maxCount){
      maxCount     = iter->second;
    }
  }

  for (iter = cMap.begin(); iter != endMap; ++iter){
    if (iter->second ==  maxCount){
      ++numberWithMax;
    }
  }

  return numberWithMax > 1;
}

short decideOnTie(const countMap& cMap) throw (std::invalid_argument) {
  short maxCount = 0;
  short numberWithMax = 0;
  int index, count;
  short temp;

  if ( cMap.empty() ) {
    // No votes to count, treat as background.
    throw std::invalid_argument("countMap should not be empty");
  }

  countMap::const_iterator iter;
  const countMap::const_iterator endMap = cMap.end();
  for (iter = cMap.begin(); iter != endMap; ++iter){
    if (iter->second > maxCount){
      maxCount     = iter->second;
    }
  }

  for (iter = cMap.begin(); iter != endMap; ++iter){
    if (iter->second ==  maxCount){
      ++numberWithMax;
    }
  }

  short tiedLabels[numberWithMax];

  count = 0;
  for (iter = cMap.begin(); iter != endMap; ++iter){
    if (iter->second ==  maxCount){
      tiedLabels[count] = iter->first;
      ++count;
    }
  }

  index = (int) floor( ran2(&ran2Seed) * count );

  temp = tiedLabels[index];
  
  return temp;
}
