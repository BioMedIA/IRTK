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

#include <irtkImageToImage2.h>

irtkImageToImage2::irtkImageToImage2()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;

  // Default parameters
  _DebugFlag = False;
}

irtkImageToImage2::~irtkImageToImage2()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;
}

void irtkImageToImage2::SetInput(irtkBaseImage*image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageToImage::SetInput: Input is not an image\n";
    exit(1);
  }
}

void irtkImageToImage2::SetOutput(irtkBaseImage *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkImageToImage::SetOutput: Output is not an image\n";
    exit(1);
  }
}

void irtkImageToImage2::Debug(const char *message)
{
  if (_DebugFlag == True) cout << message << endl;
}

void irtkImageToImage2::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageToImage::Initialize");

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }

  if (_output == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no output" << endl;
    exit(1);
  }

  if (_input->IsEmpty() == True) {
    cerr << this->NameOfClass() << "::Run: Input is empty" << endl;
    exit(1);
  }

  // Check whether filter requires buffering
  if (this->RequiresBuffering()) {
    this->Debug("irtkImageToImage::Initialize: Filter requires buffering");

    // Check whether filter has external buffer
    if (_input == _output) {
      this->Debug("irtkImageToImage::Initialize: Filter has internal buffer");
      _tmp    = _output;
      if (dynamic_cast<const irtkGenericImage<char> *>(_output) != NULL) {
        _output = new irtkGenericImage<char>(*(dynamic_cast<const irtkGenericImage<char> *>(_output)));
      } else if (dynamic_cast<const irtkGenericImage<unsigned char> *>(_output) != NULL) {
        _output = new irtkGenericImage<unsigned char>(*(dynamic_cast<const irtkGenericImage<unsigned char> *>(_output)));
      } else if (dynamic_cast<const irtkGenericImage<short> *>(_output) != NULL) {
        _output = new irtkGenericImage<short>(*(dynamic_cast<const irtkGenericImage<short> *>(_output)));
      } else if (dynamic_cast<const irtkGenericImage<unsigned short> *>(_output) != NULL) {
        _output = new irtkGenericImage<unsigned short>(*(dynamic_cast<const irtkGenericImage<unsigned short> *>(_output)));
      } else if (dynamic_cast<const irtkGenericImage<float> *>(_output) != NULL) {
        _output = new irtkGenericImage<float>(*(dynamic_cast<const irtkGenericImage<float> *>(_output)));
      } else if (dynamic_cast<const irtkGenericImage<double> *>(_output) != NULL) {
        _output = new irtkGenericImage<double>(*(dynamic_cast<const irtkGenericImage<double> *>(_output)));
      } else {
        cerr << "irtkImageToImage::Initialize: Cannot allocate image of unknown type" << endl;
        exit(1);
      }
    } else {
      this->Debug("irtkImageToImage::Initialize: Filter has external buffer");
      _tmp    = NULL;
    }
  } else {
    this->Debug("irtkImageToImage::Initialize: Filter requires no buffering");
  }

  // Make sure that output has the correct dimensions
  if (_input != _output) _output->Initialize(_input->GetImageAttributes());
}

void irtkImageToImage2::Finalize()
{
  // Print debugging information
  this->Debug("irtkImageToImage::Finalize");

  // Check whether filter requires buffering
  if (this->RequiresBuffering()) {
    this->Debug("irtkImageToImage::Finalize: Filter requires buffering");

    // Check whether filter has internal buffer
    if (_tmp != NULL) {
      this->Debug("irtkImageToImage::Finalize: Filter has internal buffer");
      // Copy buffer
      *_tmp = *_output;
      // Delete buffer
      delete _output;
      // Bend pointers back
      _output = _tmp;
    } else {
      this->Debug("irtkImageToImage::Finalize: Filter has external buffer");
    }
  } else {
    this->Debug("irtkImageToImage::Finalize: Filter requires no buffering");
  }
}

double irtkImageToImage2::Run(int, int, int, int)
{
  cerr << "Filter " << this->NameOfClass() << " has no Run(int, int, int) ";
  cerr << "member function. Using irtkImageToImage::Run." << endl;
  return 0;
}

void irtkImageToImage2::Run()
{
  int x, y, z, t;

  // Do the initial set up
  this->Initialize();

  // Calculate
  for (t = 0; t < _input->GetT(); t++) {
    for (z = 0; z < _input->GetZ(); z++) {
      for (y = 0; y < _input->GetY(); y++) {
        for (x = 0; x < _input->GetX(); x++) {
          _output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}
