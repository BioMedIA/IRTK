/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <irtkLookupTable.h>

irtkLookupTable::irtkLookupTable()
{
  lookupTable = new irtkColorRGBA[2];
  lookupTable[0] = 0;
  lookupTable[1] = 0;
  minData     = 0;
  maxData     = 1;
  minDisplay  = 0;
  maxDisplay  = 1;
  _mode = ColorMode_Luminance;
}

irtkLookupTable::~irtkLookupTable()
{
  lookupTable += minData;
  delete [] lookupTable;
}

void irtkLookupTable::Update()
{
  switch (_mode) {
  case ColorMode_Red:
    this->SetColorModeToRed();
    break;
  case ColorMode_Green:
    this->SetColorModeToGreen();
    break;
  case ColorMode_Blue:
    this->SetColorModeToBlue();
    break;
  case ColorMode_Luminance:
    this->SetColorModeToLuminance();
    break;
  case ColorMode_InverseLuminance:
    this->SetColorModeToInverseLuminance();
    break;
  case ColorMode_Rainbow:
    this->SetColorModeToRainbow();
    break;
  case ColorMode_HotMetal:
    this->SetColorModeToHotMetal();
    break;
  case ColorMode_Jacobian:
    this->SetColorModeToJacobian();
    break;
  case ColorMode_JacobianExpansion:
    this->SetColorModeToJacobianExpansion();
    break;
  case ColorMode_JacobianContraction:
    this->SetColorModeToJacobianContraction();
    break;
  case ColorMode_Custom:
    break;
  default:
    cerr << "irtkLookupTable::Update: Unknown color mode" << endl;
    exit(1);
  }
}

void irtkLookupTable::Initialize(int value1, int value2)
{
  lookupTable += minData;
  delete [] lookupTable;
  minData    = value1;
  maxData    = value2;
  minDisplay = value1;
  maxDisplay = value2;
  lookupTable  =  new irtkColorRGBA[maxData - minData + 1];
  lookupTable += -minData;
  this->Update();
}

void irtkLookupTable::SetColorModeToLuminance()
{
  int i;

  _mode = ColorMode_Luminance;
  for (i = minData; i < minDisplay; i++) {
    lookupTable[i] = 0;
    lookupTable[i].a = 1;
  }
  for (i = minDisplay; i <= maxDisplay; i++) {
    lookupTable[i] = round(((i - minDisplay) /
                            (double)(maxDisplay - minDisplay) * 255));
    lookupTable[i].a = 1;
  }
  for (i = maxDisplay+1; i < maxData; i++) {
    lookupTable[i] = 255;
    lookupTable[i].a = 1;
  }
}

void irtkLookupTable::SetColorModeToInverseLuminance()
{
  int i;

  _mode = ColorMode_InverseLuminance;
  for (i = minData; i < minDisplay; i++) {
    lookupTable[i]   = 255;
    lookupTable[i].a = 1;

  }
  for (i = minDisplay; i <= maxDisplay; i++) {
    lookupTable[i] = round(((maxDisplay - i) /
                            (double)(maxDisplay - minDisplay) * 255));
    lookupTable[i].a = 1;

  }
  for (i = maxDisplay+1; i < maxData; i++) {
    lookupTable[i]   = 0;
    lookupTable[i].a = 1;
  }
}

void irtkLookupTable::SetColorModeToHotMetal()
{
#ifdef HAS_COLOR
  int i;

  _mode = ColorMode_HotMetal;
  for (i = minData; i < minDisplay; i++) {
    lookupTable[i].r = 0;
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;

  }
  for (i = minDisplay; i <= maxDisplay; i++) {
    lookupTable[i].r = 255;
    lookupTable[i].g = round(((i - minDisplay) /
                              (double)(maxDisplay - minDisplay) * 255));
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
  for (i = maxDisplay+1; i < maxData; i++) {
    lookupTable[i].r = 255;
    lookupTable[i].g = 255;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::SetColorModeToJacobian()
{
#ifdef HAS_COLOR
  int i;

  _mode = ColorMode_Jacobian;
  for (i = minData; i < ((minDisplay < 100) ? minDisplay : 100); i++) {
    lookupTable[i].HSVtoRGB(180.0 / 360.0, 1, 1);
    lookupTable[i].a = 1;
  }
  for (i = minDisplay; i <= 100; i++) {
    lookupTable[i].HSVtoRGB(- (i - 100) / (double)(minDisplay - 100) * 60.0 / 360.0 + 240.0 / 360.0, 1, 1);
    lookupTable[i].a = 1;
  }
  for (i = 100; i <= maxDisplay; i++) {
    lookupTable[i].HSVtoRGB((i - 100) / (double)(maxDisplay - 100) * 60.0 / 360.0, 1, 1);
    lookupTable[i].a = 1;
  }
  for (i = ((maxDisplay > 100) ? maxDisplay : 100); i < maxData; i++) {
    lookupTable[i].HSVtoRGB(60.0 / 360.0, 1, 1);
    lookupTable[i].a = 1;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::SetColorModeToJacobianExpansion()
{
#ifdef HAS_COLOR
  int i, min, max;

  _mode = ColorMode_JacobianExpansion;
  min = ((minDisplay > 100) ? minDisplay : 100);
  max = ((maxDisplay > 100) ? maxDisplay : 100);
  for (i = minData; i < min; i++) {
    lookupTable[i].r = 0;
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 0;
  }
  for (i = min; i <= max; i++) {
    lookupTable[i].HSVtoRGB((i - min) / (double)(max - min + 1) * 60.0 / 360.0, 1, 1);
    lookupTable[i].a = (i - min) / (double)(max - min + 1);
  }
  for (i = max+1; i < maxData; i++) {
    lookupTable[i].HSVtoRGB(60.0 / 360.0, 1, 1);
    lookupTable[i].a = 1;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::SetColorModeToJacobianContraction()
{
#ifdef HAS_COLOR
  int i, min, max;

  _mode = ColorMode_JacobianContraction;
  min = ((minDisplay < 100) ? minDisplay : 100);
  max = ((maxDisplay < 100) ? maxDisplay : 100);
  for (i = minData; i < min; i++) {
    lookupTable[i].HSVtoRGB(180 / 360.0, 1, 1);
    lookupTable[i].a = 1;
  }
  for (i = min; i <= max; i++) {
    lookupTable[i].HSVtoRGB((i - max) / (double)(max - min + 1) * 60.0 / 360.0 + 240.0 / 360.0, 1, 1);
    lookupTable[i].a = -(i - max) / (double)(max - min + 1);
  }
  for (i = max+1; i < maxData; i++) {
    lookupTable[i].r = 0;
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 0;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::SetColorModeToRed()
{
#ifdef HAS_COLOR
  int i;

  _mode = ColorMode_Red;
  for (i = minData; i < minDisplay; i++) {
    lookupTable[i].r = 0;
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
  for (i = minDisplay; i <= maxDisplay; i++) {
    lookupTable[i].r = round(((i - minDisplay) /
                              (double)(maxDisplay - minDisplay) * 255));
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
  for (i = maxDisplay+1; i < maxData; i++) {
    lookupTable[i].r = 255;
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::SetColorModeToGreen()
{
#ifdef HAS_COLOR
  int i;

  _mode = ColorMode_Green;
  for (i = minData; i < minDisplay; i++) {
    lookupTable[i].r = 0;
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
  for (i = minDisplay; i <= maxDisplay; i++) {
    lookupTable[i].g = round(((i - minDisplay) /
                              (double)(maxDisplay - minDisplay) * 255));
    lookupTable[i].r = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
  for (i = maxDisplay+1; i < maxData; i++) {
    lookupTable[i].g = 255;
    lookupTable[i].r = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::SetColorModeToBlue()
{
#ifdef HAS_COLOR
  int i;

  _mode = ColorMode_Blue;
  for (i = minData; i < minDisplay; i++) {
    lookupTable[i].r = 0;
    lookupTable[i].g = 0;
    lookupTable[i].b = 0;
    lookupTable[i].a = 1;
  }
  for (i = minDisplay; i <= maxDisplay; i++) {
    lookupTable[i].b = round(((i - minDisplay) /
                              (double)(maxDisplay - minDisplay) * 255));
    lookupTable[i].g = 0;
    lookupTable[i].r = 0;
    lookupTable[i].a = 1;
  }
  for (i = maxDisplay+1; i < maxData; i++) {
    lookupTable[i].b = 255;
    lookupTable[i].g = 0;
    lookupTable[i].r = 0;
    lookupTable[i].a = 1;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::SetColorModeToRainbow()
{
#ifdef HAS_COLOR
  int i;

  _mode = ColorMode_Rainbow;
  for (i = minData; i < minDisplay; i++) {
    lookupTable[i].HSVtoRGB(2.0/3.0, 1, 1);
    lookupTable[i].a = 1;
  }
  for (i = minDisplay; i <= maxDisplay; i++) {
    lookupTable[i].HSVtoRGB((maxDisplay - i) / (double)(maxDisplay - minDisplay) * 2.0/3.0, 1, 1);
    lookupTable[i].a = 1;
  }
  for (i = maxDisplay+1; i < maxData; i++) {
    lookupTable[i].HSVtoRGB(0, 1, 1);
    lookupTable[i].a = 1;
  }
#else
  this->SetColorModeToLuminance();
#endif
}

void irtkLookupTable::Read(char *filename)
{
  float a;
  char buffer[255];
  int i, r, g, b, value1, value2;

  // Open file
  ifstream from(filename);

  // Check file format
  from >> buffer;
  if (strcmp(buffer, "irtkLookupTable") != 0) {
    cerr << "irtkLookupTable::Read: Can't read lookup table" << endl;
    exit(1);
  }

  // Read min and max values
  from >> value1 >> value2;

  // Initialize lookup table
  this->Initialize(value1, value2);

  for (i = minData; i <= maxData; i++) {
    from >> r >> g >> b >> a;
    lookupTable[i].r = r;
    lookupTable[i].g = g;
    lookupTable[i].b = b;
    lookupTable[i].a = a;
  }

  // This is a custom lookup table
  _mode = ColorMode_Custom;
}

void irtkLookupTable::Write(char *filename)
{
  cerr << "irtkLookupTable::Write: Not yet implemented" << endl;
  exit(1);
}

