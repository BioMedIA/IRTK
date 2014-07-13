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

#ifndef combineLabels_core_h
#define combineLabels_core_h

#include <map>
#include <stdexcept>

#include <boost/random/mersenne_twister.hpp>

typedef std::map <short, short> countMap;

short int getMostPopular (const countMap& cMap);
bool 	  isEquivocal 	 (const countMap& cMap);
short int decideOnTie 	 (const countMap& cMap, boost::mt19937& rng) throw (std::invalid_argument);

#endif // combineLabels_core_h
