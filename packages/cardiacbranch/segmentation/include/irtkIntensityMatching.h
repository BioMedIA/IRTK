/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKINTENSITYMATCHING_H

#define _IRTKINTENSITYMATCHING_H

#include <irtkGeometry.h>

#include <irtkImage.h>

class irtkIntensityMatching : public irtkObject
{

protected:

///minimum parameter
  double _min;
///maximum parameter
  double _max;
///minimum parameter
  double _ymin;
///maximum parameter
  double _ymax;
///best fit line
  double _m;
///best git line
  double _b;

  irtkGreyImage _hist;
  irtkGreyImage _trans;

public:

  /// Constructor
  irtkIntensityMatching();

  /// Constructor
  irtkIntensityMatching(double, double, double, double);

  void SetMinMax(double min, double max);
  void SetYMinMax(double min, double max);
  void Proj(double& x, double& y);

  /// Copy Constructor
  irtkIntensityMatching(const irtkIntensityMatching &);

  /// Destructor
  ~irtkIntensityMatching();

  void PCA(double *x, double *y, double *weights, int no);
  void WeightedPCA(double *x, double *y, double *weights, int no);
  void WeightedLeastSqaures(double *x1, double *y1, double *weights, int no);
  void MatchMeanAndVariance(double *x1, double *y1, double *weights, int no);

  double Lin(double x);


  /// Print info
  void Print();
  void NameOfClass();

  void WriteHist(char *);
  void WriteTrans(char *);
};

inline void irtkIntensityMatching::SetMinMax(double min, double max)
{
  _min=min;
  _max=max;
}

inline void irtkIntensityMatching::SetYMinMax(double min, double max)
{
  _ymin=min;
  _ymax=max;
}

inline void irtkIntensityMatching::WriteHist(char *name)
{
  _hist.Write(name);
}

inline void irtkIntensityMatching::WriteTrans(char *name)
{
  _trans.Write(name);
}

inline void irtkIntensityMatching::NameOfClass()
{
  cerr<<"irtkIntensityMatching"<<endl;
}

inline void irtkIntensityMatching::Proj(double& x, double& y)
{
  double norm=sqrt(1+_m*_m);
  double xx,yy;

  xx=(x+_m*y)/norm;
  yy=(_m*x-y+_b)/norm;
  x=xx;
  y=yy;
}

inline double irtkIntensityMatching::Lin(double x)
{
  return _m*x+_b;
}



#endif
