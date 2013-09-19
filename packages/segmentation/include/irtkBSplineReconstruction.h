/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstruction.h 147 2010-05-03 17:41:54Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2010-05-03 18:41:54 +0100 (Mon, 03 May 2010) $
  Version   : $Revision: 147 $
  Changes   : $Author: mm3 $

=========================================================================*/

#ifndef _irtkBSplineReconstruction_H

#define _irtkBSplineReconstruction_H

#define LOOKUPTABLESIZE 1000
#define DBL_LUTSIZE (double)(LOOKUPTABLESIZE-1)
#define MAX_IMAGES 500


#include <irtkImage.h>
#include <irtkResampling.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

#include <vector>
using namespace std;


class irtkBSplineReconstruction : public irtkObject
{

protected:
  
irtkRealImage **_dx;
irtkRealImage **_coeffs;
irtkRealImage **_weights;

irtkRealImage *_temp;
irtkGreyImage *_temp2;

irtkRealImage *_template;
irtkRealImage *_evalTarget;

//irtkRealImage **_input      = NULL;
//irtkTransformation **_transf = NULL;

int _inputCount;

double _padding;

double LookupTable  [LOOKUPTABLESIZE][4];

inline double Bsp(int i, double t);
void initialiseLookupTable();
void clearRealImage(irtkRealImage *img);
void clearGreyImage(irtkGreyImage *img);
irtkRealImage * zeroPad(irtkRealImage *in, int count);
irtkRealImage * growImage(irtkRealImage *in, int levels);
irtkRealImage * downsampleFactorTwo(irtkRealImage *in);
void subdivide(irtkRealImage *in, irtkRealImage *out);
double getIntensity(irtkRealImage *coeffs, double xin, double yin, double zin);
void evaluateCoeffsOnImageLattice(irtkRealImage *img, irtkRealImage *coeffs);
void evaluateCoeffsOnImageLattice(irtkGreyImage *img, irtkRealImage *coeffs);
void estimateCoeffs(int res, vector<irtkRealImage>& _input, vector<irtkRigidTransformation>& _transf);
void conv_3D(irtkRealImage *image, irtkRealImage *temp, double *ker, int lKern);
void conv_3D_long(irtkRealImage *image, irtkRealImage *temp, double *ker, int lKern);
double getIntensity_verbose(irtkRealImage *coeffs, double x, double y, double z);
irtkRealImage * upsampleFactorTwo(irtkRealImage *in);


public:
  
  irtkBSplineReconstruction();
  
  void Reconstruct(int levels, int lastLevel, irtkRealImage& temp, vector<irtkRealImage>& _input, vector<irtkRigidTransformation>& _transf);

};

inline double irtkBSplineReconstruction::Bsp(int i, double t)
{
  switch (i) {
  case 0:
    return (1-t)*(1-t)*(1-t)/6.0;
  case 1:
    return (3*t*t*t - 6*t*t + 4)/6.0;
  case 2:
    return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
  case 3:
    return (t*t*t)/6.0;
  }
  return 0;
}

#endif
