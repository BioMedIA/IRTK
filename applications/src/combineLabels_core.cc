/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

// avoid tweaking build system for the moment
#include "../include/combineLabels_core.h"
#include <irtkImage.h>

#include <sys/types.h>
#include <stdexcept>

#ifdef WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
              

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

short decideOnTie(const countMap& cMap, boost::random::mt19937& rng) throw (std::invalid_argument) {
  short maxCount = 0;
  short numberWithMax = 0;
  int   count;

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

  boost::uniform_int<> uniformCount(0, count -1);
  boost::variate_generator<boost::mt19937&, boost::uniform_int<> > indexEngine(rng, uniformCount);

  return tiedLabels[indexEngine()];
}
