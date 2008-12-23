/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKLOOKUPTABLE_H

#define _IRTKLOOKUPTABLE_H

#include <irtkColorRGBA.h>
#include <irtkColor.h>

typedef enum { ColorMode_Custom,
               ColorMode_HotMetal,
               ColorMode_Red,
               ColorMode_Green,
               ColorMode_Blue,
               ColorMode_Jacobian,
               ColorMode_JacobianExpansion,
               ColorMode_JacobianContraction,
               ColorMode_Luminance,
               ColorMode_InverseLuminance,
               ColorMode_Rainbow
             } irtkColorMode;


class irtkLookupTable
{

public:

  /// Color mode
  irtkColorMode _mode;

  /// Min and max value of data
  int minData;
  int maxData;

  /// Min and max value of display
  int minDisplay;
  int maxDisplay;

  /// Lookup table
  irtkColorRGBA *lookupTable;

  /// Constructor
  irtkLookupTable();

  /// Destructor
  ~irtkLookupTable();

  /// Update lookup table
  void Update();

  /// Initialize lookup table with min and max values
  void Initialize(int, int);

  /// Set minimum and maximum intensities
  void SetMinMaxIntensity(int,   int);

  /// Get minimum and maximum intensities
  void GetMinMaxIntensity(int &, int &);

  /// Set minimum intensity
  void SetMinIntensity(int);

  /// Get minimum intensity
  int  GetMinIntensity();

  /// Set maximum intensity
  void SetMaxIntensity(int);

  /// Get maximum intensity
  int  GetMaxIntensity();

  /// Color scheme functions
  void SetColorModeToLuminance();
  void SetColorModeToInverseLuminance();
  void SetColorModeToRed();
  void SetColorModeToGreen();
  void SetColorModeToBlue();
  void SetColorModeToRainbow();
  void SetColorModeToHotMetal();
  void SetColorModeToJacobian();
  void SetColorModeToJacobianExpansion();
  void SetColorModeToJacobianContraction();

  /// Return color scheme
  irtkColorMode GetColorMode();

  /// Read lookup table from file
  void Read(char *);

  /// Write lookup table to file
  void Write(char *);

};

inline int irtkLookupTable::GetMinIntensity()
{
  return minDisplay;
}

inline int irtkLookupTable::GetMaxIntensity()
{
  return maxDisplay;
}

inline void irtkLookupTable::SetMinIntensity(int value)
{
  if (value > maxData) {
    minDisplay = maxData;
  } else {
    minDisplay = value;
  }
  this->Update();
}

inline void irtkLookupTable::SetMaxIntensity(int value)
{
  if (value < minData) {
    maxDisplay = minData;
  } else {
    maxDisplay = value;
  }
  this->Update();
}

inline void irtkLookupTable::SetMinMaxIntensity(int value1, int value2)
{
  if (value1 > maxData) {
    minDisplay = maxData;
  } else {
    minDisplay = value1;
  }
  this->Update();
  if (value2 < minData) {
    maxDisplay = minData;
  } else {
    maxDisplay = value2;
  }
  this->Update();
}

inline irtkColorMode irtkLookupTable::GetColorMode()
{
  return _mode;
}


#endif
