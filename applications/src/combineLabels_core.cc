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


char *output_name = NULL, **input_names = NULL;
char *mask_name = NULL;

long ran2Seed;
long ran2initialSeed;

// first short is the label, second short is the number of images voting
// for that label.
// typedef map<short, short> countMap;

map<short, short>::iterator iter;

short getMostPopular(countMap cMap){
  short maxCount = 0, mostPopLabel = -1;

  if (cMap.size() == 0){
    // No votes to count, treat as background.
    return 0;
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){

    if (iter->second > maxCount){
      maxCount     = iter->second;
      mostPopLabel = iter->first;
    }
  }

  return mostPopLabel;
}

bool isEquivocal(countMap cMap){
  short maxCount = 0;
  short numberWithMax = 0;

  if (cMap.size() == 0){
    // No votes to count, treat as background.
    return false;
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second > maxCount){
      maxCount     = iter->second;
    }
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second ==  maxCount){
      ++numberWithMax;
    }
  }

  if (numberWithMax > 1){
    return true;
  } else {
    return false;
  }

}

short decideOnTie(countMap cMap) throw (std::invalid_argument) {
  short maxCount = 0;
  short numberWithMax = 0;
  int index, count;
  short temp;

  if (cMap.empty()){
    // No votes to count, treat as background.
    throw std::invalid_argument("countMap should not be empty");
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second > maxCount){
      maxCount     = iter->second;
    }
  }

  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second ==  maxCount){
      ++numberWithMax;
    }
  }

  short *tiedLabels = new short[numberWithMax];

  count = 0;
  for (iter = cMap.begin(); iter != cMap.end(); ++iter){
    if (iter->second ==  maxCount){
      tiedLabels[count] = iter->first;
      ++count;
    }
  }

  index = (int) floor( ran2(&ran2Seed) * count );

  temp = tiedLabels[index];
  delete tiedLabels;

  return temp;
}
