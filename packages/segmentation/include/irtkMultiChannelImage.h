/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMULTICHANNELIMAGE_H

#define _IRTKMULTICHANNELIMAGE_H

#include <irtkObject.h>
#include <irtkImage.h>
#include <vector>

using namespace std;

/**

Atlas probability map class

*/

class irtkMultiChannelImage : public irtkObject
{

protected:

  //Vector of all images
  vector<irtkRealImage> _images;
  //Vector of pointers - current position in images
  vector<irtkRealPixel*> _pointers;

  ///Number of voxels
  int _number_of_voxels;
  int _padding;

  irtkRealImage _mask;

public:

  ///Adds new image (channel)
  void AddImage(const irtkRealImage &image);

  ///Returns reference to image (channel)
  irtkRealImage& GetImage(int channel);

  void SetImage(int channel, irtkRealImage &image);

  /// Moves pointers in all images to the first voxel
  void First();

  /// Moves pointers in all images to the next voxel
  void Next();

  ///Returns intensity value at pointer
  irtkRealPixel GetValue(int channel);

  ///Returns intensity value at position
  irtkRealPixel GetValue(int x, int y, int z, int channel);

  ///Sets intensity value at pointer
  void SetValue(int tissue, irtkRealPixel value);

  ///Sets intensity value at position
  void SetValue(int x, int y, int z, int tissue, irtkRealPixel value);

  /// Returns number of voxels
  int GetNumberOfVoxels();
  int GetX();
  int GetY();
  int GetZ();

  ///Saves image to a file
  void Write(int channel, const char *file_name);

  ///Returns number of channels
  int GetNumberOfChannels();

  ///Returns minimum and maximum intenzity of channel

  void GetMinMax(int channel, int &min, int &max);

  ///Retuns average of all images
  irtkRealImage Average();
  ///Retuns sum of all images
  irtkRealImage Add();
  /// Returns the product of all images
  irtkRealImage Multiply();
  double ComputeVolume(int channel, int label = 1);
  irtkRealImage Subtract();
  irtkRealImage Divide(double scale);
  void Log(int channel, double scale=1000);
  void Exp(int channel, double scale=1000);
  irtkRealImage Max();
  irtkRealImage CreateMask();
  void Brainmask();
  void SetPadding(int);
  void HistogramMatching(int tchannel, int rchannel, double sigma=0);
  void HistogramEqualization(int channel);
  void AdjustMean(int tchannel, int rchannel);
};

#endif
