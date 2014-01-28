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

#include <irtkImage.h>

#include <irtkImageFunction.h>

irtkImageFunction::irtkImageFunction()
{
  // Set input
  _input  = NULL;

  // Default parameters
  _DebugFlag    = false;
  _DefaultValue = 0;
}

irtkImageFunction::~irtkImageFunction()
{
  // Set input
  _input  = NULL;
}

void irtkImageFunction::SetInput(irtkImage *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageFunction::SetInput: Input is not an image\n";
    exit(1);
  }
}

void irtkImageFunction::Debug(const char *message)
{
  if (_DebugFlag == true) cout << message << endl;
}

void irtkImageFunction::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageFunction::Initialize");

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }
}

