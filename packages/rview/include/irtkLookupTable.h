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

	friend class irtkRView;
	
  /// Min value of data
  int _minData;

  /// Max value of data
  int _maxData;

  /// Min value of display
  int _minDisplay;

  /// Max value of display
  int _maxDisplay;

  /// Color mode
  irtkColorMode _mode;

  /// Update lookup table
  void Update();

  /// Initialize lookup table with min and max values
  void Initialize(int, int);

  /// Set minimum and maximum display intensities
  void SetMinMaxDisplayIntensity(int,   int);

  /// Get minimum and maximum display intensities
  void GetMinMaxDisplayIntensity(int &, int &);

  /// Set minimum display intensity
  void SetMinDisplayIntensity(int);

  /// Set maximum display intensity
  void SetMaxDisplayIntensity(int);

  /// Get maximum display intensity
  int  GetMaxDisplayIntensity();

  /// Get minimum intensity
  int  GetMinIntensity();

  /// Get maximum intensity
  int  GetMaxIntensity();

 
public:

  /// Lookup table
  irtkColorRGBA *lookupTable;

  /// Constructor
  irtkLookupTable(int = 0, int = 1);

  /// Destructor
  ~irtkLookupTable();

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

  /// Get minimum display intensity
  int  GetMinDisplayIntensity();

  /// Return color scheme
  irtkColorMode GetColorMode();

  /// Read lookup table from file
  void Read(char *);

  /// Write lookup table to file
  void Write(char *);

};

inline int irtkLookupTable::GetMinIntensity()
{
  return _minData;
}

inline int irtkLookupTable::GetMaxIntensity()
{
  return _maxData;
}

inline int irtkLookupTable::GetMinDisplayIntensity()
{
  return _minDisplay;
}

inline int irtkLookupTable::GetMaxDisplayIntensity()
{
  return _maxDisplay;
}

inline void irtkLookupTable::SetMinDisplayIntensity(int value)
{
  if (value > _maxData) {
    _minDisplay = _maxData;
  } else {
    _minDisplay = value;
  }
  this->Update();
}

inline void irtkLookupTable::SetMaxDisplayIntensity(int value)
{
  if (value < _minData) {
    _maxDisplay = _minData;
  } else {
    _maxDisplay = value;
  }
  this->Update();
}

inline void irtkLookupTable::SetMinMaxDisplayIntensity(int value1, int value2)
{
  if (value1 > _maxData) {
    _minDisplay = _maxData;
  } else {
    _minDisplay = value1;
  }
  this->Update();
  if (value2 < _minData) {
    _maxDisplay = _minData;
  } else {
    _maxDisplay = value2;
  }
  this->Update();
}

inline irtkColorMode irtkLookupTable::GetColorMode()
{
  return _mode;
}


#endif
