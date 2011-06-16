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

#include <irtkProbabilisticAtlas.h>

irtkProbabilisticAtlas::irtkProbabilisticAtlas()
{
  _number_of_voxels=0;
  _number_of_tissues=0;
}

void irtkProbabilisticAtlas::AddImage(irtkRealImage image)
{
  if (_images.size() == 0) {
    _number_of_voxels = image.GetNumberOfVoxels();
  } else {
    if (_number_of_voxels != image.GetNumberOfVoxels()) {
      cerr << "Image sizes mismatch" << endl;
      exit(1);
    }
  }
  _images.push_back(image);
  _pointers.push_back(image.GetPointerToVoxels());
  _number_of_tissues = _images.size();
}

void irtkProbabilisticAtlas::NormalizeAtlas(irtkRealImage background)
{
  int i, j;
  irtkRealPixel norm;

  // Add extra image
  if (_number_of_tissues == 0) {
    cerr << "irtkProbabilisticAtlas::NormalizeAtlas: No probability maps found" << endl;
    exit(1);
  } else {
    this->AddImage(background);
    _number_of_tissues = _images.size();

  }

  // normalize atlas to 0 to 1
  this->First();
  for (i = 0; i < _number_of_voxels; i++) {
    norm = 0;
    for (j = 0; j < _number_of_tissues; j++) {
      norm += *(_pointers[j]);
    }
    if (norm>0) {
      for (j = 0; j < _number_of_tissues; j++)
        *(_pointers[j]) = *(_pointers[j])/norm;
    } else {
      for (j = 0; j < _number_of_tissues-1; j++)
        *(_pointers[j]) = 0;
      *(_pointers[_number_of_tissues-1]) = 1;
    }
    this->Next();
  }
}

void irtkProbabilisticAtlas::NormalizeAtlas()
{
  int i, j;
  irtkRealPixel norm;

  // Add extra image
  if (_number_of_tissues == 0) {
    cerr << "irtkProbabilisticAtlas::NormalizeAtlas: No probability maps found" << endl;
    exit(1);
  }
  _number_of_tissues = _images.size();


  // normalize atlas to 0 to 1
  this->First();
  for (i = 0; i < _number_of_voxels; i++) {
    norm = 0;
    for (j = 0; j < _number_of_tissues; j++) {
      norm += *(_pointers[j]);
    }
    if (norm>0) {
      for (j = 0; j < _number_of_tissues; j++)
        *(_pointers[j]) = *(_pointers[j])/norm;
    } else {
      for (j = 0; j < _number_of_tissues-1; j++)
        *(_pointers[j]) = 0;
      *(_pointers[_number_of_tissues-1]) = 1;
    }
    this->Next();
  }
}


void irtkProbabilisticAtlas::AddBackground()
{
  int i, j;
  irtkRealPixel norm, min, max;

  // Add extra image
  if (_number_of_tissues == 0) {
    cerr << "irtkProbabilisticAtlas::AddBackground: No probability maps found" << endl;
    exit(1);
  } else {
    this->AddImage(_images[0]);
    _number_of_tissues = _images.size();

  }

  // normalize atlas to 0 to 1
  irtkRealImage other;
  other = _images[0];
  for (i = 1; i < _number_of_tissues-1; i++) {
    cerr<<i<<endl;
    other += _images[i];
  }
  other.GetMinMax(&min, &max);

  this->First();
  for (i = 0; i < _number_of_voxels; i++) {
    norm = 0;
    for (j = 0; j < _number_of_tissues-1; j++) {
      *(_pointers[j]) = (*(_pointers[j]) - min) / (max - min);
      norm += *(_pointers[j]);
    }
    double value = 1 - norm;
    if (value < 0) value = 0;
    if (value > 1) value = 1;
    *(_pointers[_number_of_tissues-1]) = value;
    this->Next();
  }
}

void irtkProbabilisticAtlas::AddProbabilityMaps(int n, irtkRealImage **atlas)
{
  int i;

  for (i = 0; i < n; i++) {
    this->AddImage(*(atlas[i]));
  }
};

void irtkProbabilisticAtlas::Write(int i, const char *filename)
{
  if  (i < _number_of_tissues) {
    _images[i] *= 255;
    _images[i].Write(filename);
    _images[i] /= 255;
  } else {
    cerr << "irtkProbabilisticAtlas::Write: No such probability map" << endl;
    exit(1);
  }
}

irtkRealImage irtkProbabilisticAtlas::ComputeHardSegmentation()
{
  cerr<<"irtkProbabilisticAtlas::ComputeHardSegmentation"<<endl;
  int i, tissue, j=0;
  double max = 0;
  _segmentation = _images[0];
  First();
  irtkRealPixel *ptr = _segmentation.GetPointerToVoxels();

  for (i = 0; i < _number_of_voxels; i++) {
    max=0;
    tissue = -1;
    for (j=0; j< GetNumberOfTissues(); j++) {
      if (GetValue(j) > max) {
        tissue = j;
        max = GetValue(j);
      }
    }

    /*    switch (tissue)
        {
          case 0: *ptr = -1; break;
          case 1: *ptr = 45; break;
          case 2: *ptr = 75; break;
          case 3: *ptr = 80; break;
          case 4: *ptr = 82; break;
          case 5: *ptr = 85; break;
          case 6: *ptr = 81; break;
          case 7: *ptr = 87; break;
          case 8: *ptr = 89; break;
          case 9: *ptr = 94; break;
          case 10: *ptr = 100; break;
          case 11: *ptr = 100; break;

          default: cerr << "Tissue identificator out of range." <<endl;
        }
    */
    *ptr = tissue;
    Next();
    ptr++;

    /*switch (tissue)
    {
      case WHITE: *ptr = white_value; break;
      case GREY:  *ptr = grey_value;  break;
      case CSF:   *ptr = csf_value;   break;
      case OTHER: *ptr = 0; break;

      default: cerr << "Tissue identificator out of range." <<endl;
    */
  }
  return _segmentation;
}

irtkRealImage irtkProbabilisticAtlas::GetImage(int i)
{
  return _images[i];
}

void irtkProbabilisticAtlas::ExtractLabel(int label, irtkRealImage& image)
{
  cerr<<"irtkProbabilisticAtlas::ExtractLabel"<<endl;
  int i;
  //ComputeHardSegmentation();

  irtkRealPixel *ptr = _segmentation.GetPointerToVoxels();
  irtkRealPixel *ptr_im = image.GetPointerToVoxels();

  for (i = 0; i < _number_of_voxels; i++) {
    if (*ptr == label) *ptr_im = 1;
    else *ptr_im = 0;
    ptr++;
    ptr_im++;
  }
}

void irtkProbabilisticAtlas::WriteHardSegmentation(const char *filename)
{
  _segmentation.Write(filename);
}
