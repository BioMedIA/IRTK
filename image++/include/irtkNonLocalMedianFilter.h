/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKNONLOCALMEDIANFILTER_H

#define _IRTKNONLOCALMEDIANFILTER_H

#include <irtkImageToImage.h>

/**
 * Class for non local median filter of images
 *
 * This class defines and implements the non local median filter of images. 
 * Follwing Y. Li and S. Osher. A new median formula with applications to PDE
 * based denoising. Commun. Math. Sci., 7(3):741–753, 2009.
 * The non local weight is a Gaussian function of 4D distance where intensity is
 * considered as the fourth dimension
 */

template <class VoxelType> class irtkNonLocalMedianFilter : public irtkImageToImage<VoxelType>
{

protected:

  /// Sigma (standard deviation of Gaussian kernel)
  int _Sigma;

  /// Second input, i.e. the scale image
  irtkGenericImage<irtkGreyPixel> *_input2;

  /// Third input, i.e. occlusion
  irtkGenericImage<irtkRealPixel> *_input3;

  /// Fourth input, i.e. constrain
  irtkGenericImage<VoxelType> *_input4;

  irtkGenericImage<VoxelType> *_edge;

  float *_localweight;

  float *_localneighbor;

  float _dx;

  float _dy;

  float _dz;

  float _ds;

  float _Lambda;

  /// Initialize the filter
  virtual void Initialize();

  /// Finalize the filter
  virtual void Finalize();

  virtual double EvaluateWeight(const double &);

  /** Run This method is protected and should only
  *  be called from within public member function Run().
  */
  virtual double Run(int, int, int, int);

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor to remvoe the non-local term
  /// only supply the size argument
  /// otherwise an non-local weight will be cauculated using input1
  /// i.e. second argument
  irtkNonLocalMedianFilter(int,irtkGenericImage<irtkGreyPixel>* = NULL,
      irtkGenericImage<irtkRealPixel>* = NULL,irtkGenericImage<VoxelType>* = NULL);

  /// Destructor
  ~irtkNonLocalMedianFilter();

  /// Run
  virtual void Run();

  /// Set displacement
  SetMacro(input2, irtkGenericImage<irtkGreyPixel>*);

  /// Set occlusion
  SetMacro(input3, irtkGenericImage<irtkRealPixel>*);

  /// Set constrain
  SetMacro(input4, irtkGenericImage<VoxelType>*);

  /// Set sigma
  SetMacro(Sigma, int);

  /// Set sigma
  SetMacro(Lambda, float);

};

template <class VoxelType> inline double irtkNonLocalMedianFilter<VoxelType>::EvaluateWeight(const double &distancev)
{
    return exp(- distancev / (_Sigma*_Sigma/2.0));
}

#endif
