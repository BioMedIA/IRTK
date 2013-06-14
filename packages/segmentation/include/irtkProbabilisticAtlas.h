/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPROBABILISTICATLAS_H

#define _IRTKPROBABILISTICATLAS_H

#include <irtkObject.h>

#include <irtkImage.h>

#include <vector>


/**

Atlas probability map class

*/

class irtkProbabilisticAtlas : public irtkObject
{
  // Vector of probability maps
  vector<irtkRealImage> _images;

  // Vector of pointers to probability maps
  vector<irtkRealPixel *> _pointers;

  ///Number of voxels
  int _number_of_voxels;

  ///Number of voxels
  int _number_of_tissues;

  /// Hard segmentation
  irtkRealImage _segmentation;


public:

  /// Constructor
  irtkProbabilisticAtlas();

  /// replaces image a with image
  void ReplaceImage(int, irtkRealImage);

  /// swaps images within the prob atlas
  void SwapImages(int, int);

  /// Adds new image (channel)
  void AddImage(irtkRealImage image);

  /// Moves pointers in all images to the first voxel
  void First();

  /// Moves pointers in all images to the next voxel
  void Next();

  ///Returns intensity value at pointer
  irtkRealPixel GetValue(unsigned int channel);

  ///Returns intensity value at position
  irtkRealPixel GetValue(int x, int y, int z, unsigned int tissue);

  /// Returns intensity value at position
  irtkRealPixel GetValue(int x, int y, int z, int t, unsigned int tissue);

  ///Sets intensity value at pointer
  void SetValue(unsigned int tissue, irtkRealPixel value);

  ///Sets intensity value at position
  void SetValue(int x, int y, int z, unsigned int tissue, irtkRealPixel value);

  /// Sets intensity value at position
  void SetValue(int x, int y, int z, int t, unsigned int tissue, irtkRealPixel value);

  /// Returns number of voxels
  int GetNumberOfVoxels();

  /// Add probability maps
  void AddProbabilityMaps(int, irtkRealImage **atlas);

  /// Computes probability map for the last tissue type
  void AddBackground();

  /// normalize atlas without adding background
  void NormalizeAtlas();

  /// normalize atlas while adding background
  void NormalizeAtlas(irtkRealImage background);

  ///Returns number of tissues
  int GetNumberOfTissues();

  /// Write
  void Write(int, const char *);

  /// Computes hard segmentation
  irtkRealImage ComputeHardSegmentation();

  /// Extract a label from hard segmentation
  void ExtractLabel(int label, irtkRealImage& image);

  /// Writes hard segmentation into a file
  void WriteHardSegmentation(const char *filename);

  irtkRealImage GetImage(int);

};

inline void irtkProbabilisticAtlas::First()
{
  unsigned int i;
  for (i=0; i<_pointers.size(); i++) _pointers[i] = _images[i].GetPointerToVoxels();
}

inline void irtkProbabilisticAtlas::Next()
{
  unsigned int i;
  for (i=0; i<_pointers.size(); i++) _pointers[i]++;
}

inline irtkRealPixel irtkProbabilisticAtlas::GetValue(unsigned int channel)
{
  if (channel < _images.size()) return *_pointers[channel];
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}

inline irtkRealPixel irtkProbabilisticAtlas::GetValue(int x, int y, int z, unsigned int channel)
{
  if (channel < _images.size()) return _images[channel].Get(x,y,z);
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}

inline irtkRealPixel irtkProbabilisticAtlas::GetValue(int x, int y, int z, int t, unsigned int channel)
{
  if (channel < _images.size()) return _images[channel].Get(x,y,z,t);
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}

inline void irtkProbabilisticAtlas::SetValue(unsigned int channel, irtkRealPixel value)
{
  if (channel < _images.size()) *_pointers[channel] = value;
  else {
    cerr << "Channel identificator " << channel << " out of range." <<endl;
    exit(1);
  }
}

inline void irtkProbabilisticAtlas::SetValue(int x, int y, int z, unsigned int channel, irtkRealPixel value)
{
  if (channel < _images.size()) _images[channel].Put( x, y, z, value);
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}

inline void irtkProbabilisticAtlas::SetValue(int x, int y, int z, int t, unsigned int channel, irtkRealPixel value)
{
  if (channel < _images.size()) _images[channel].Put( x, y, z, t, value);
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}

inline int irtkProbabilisticAtlas::GetNumberOfVoxels()
{
  return _number_of_voxels;
}

inline int irtkProbabilisticAtlas::GetNumberOfTissues()
{
  return _number_of_tissues;
}

#endif
