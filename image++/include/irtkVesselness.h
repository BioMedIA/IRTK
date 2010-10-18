#ifndef _IRTKVESSEL_H

#define _IRTKVESSEL_H

//#include <irtkHessianImage.h>

template <class VoxelType> class irtkVessel : public irtkObject
{

private:

  /// Debugging flag
  bool _DebugFlag;
  double _alpha;
  double _beta;
  double _c;
  double _r;


  double _sigmaMin;
  double _sigmaMax;
  double _vesselness;
  static double _sigma;
  irtkVector reference;

protected:

  /// Input image
  static irtkGenericImage<VoxelType>* _input;

  irtkGenericImage<VoxelType> _imageXX;
  irtkGenericImage<VoxelType> _imageXY;
  irtkGenericImage<VoxelType> _imageXZ;
  irtkGenericImage<VoxelType> _imageYY;
  irtkGenericImage<VoxelType> _imageYZ;
  irtkGenericImage<VoxelType> _imageZZ;


public:
  static irtkGenericImage<VoxelType> _probability;
  irtkHessianImage<VoxelType> hessianimage;
  static irtkGenericImage<VoxelType> _sigmaImage;
  irtkRealImage _vesselnessImage;
  static irtkGenericImage<VoxelType> _hxx;
  static irtkGenericImage<VoxelType> _hxy;
  static irtkGenericImage<VoxelType> _hxz;
  static irtkGenericImage<VoxelType> _hyy;
  static irtkGenericImage<VoxelType> _hyz;
  static irtkGenericImage<VoxelType> _hzz;
  static irtkGenericImage<VoxelType> _dx;
  static irtkGenericImage<VoxelType> _dy;
  static irtkGenericImage<VoxelType> _dz;

  /// Constructor
  irtkVessel();


  /// Deconstuctor
  virtual ~irtkVessel();

  void SetReference(double, double, double);
  void ChangeReference();

  /// Set input image
  void SetInput (irtkGenericImage<VoxelType> *);

  void SetAlpha (double);

  double GetAlpha();

  void SetBeta (double);

  double GetBeta();

  void SetC (double);

  double GetC();

  void SetR(double);

  double GetR();

  void SetSigma (double);

  static double GetSigma(double, double, double);

  void SetSigmaRange(double, double);


  void Initialize();

  static void increasingMagnitudeSort(int n, irtkVector& eigval, irtkMatrix& eigvec);

  void GetHessian(irtkMatrix, double, double, double);

  ///Evaluate the vesselness image
  double Vesselness(double, double, double);

  //  double ComputeVesselness(irtkVector&, double, double, double);

  void Vesselness();
  void VesselnessFDM();
  //double GetSigma(double x, double y, double z);

  void LorenzFilter();

  void SatoFilter();
  static float Ridgeness(float p[3]);
  static float Ridgeness2(float p[3]);
  static void Dridgeness(float p[3], float derivative[3]);
  static void Dridgeness2(float p[3], float derivative[3]);


  void Traversal(float p[3]);
  void Traversal2(float p[3]);


  //  void Extract(float*,float*, float*, float stepsize);
  void Extract(float*, float*, float stepsize);

  //gaussian error function
  double erf(double);

  //complementary error function
  double erfc(double);
  void IntensityMeasure(irtkGenericImage<VoxelType> *,irtkGenericImage<VoxelType> *, double, double, double);

  void CostImage(irtkGenericImage<VoxelType> *, irtkGenericImage<VoxelType> *, irtkGenericImage<VoxelType> *, double epsilon);


  /// Returns the name of the class
  char *NameOfClass();

  /// Set debugging flag
  SetMacro(DebugFlag, bool);

  /// Get debugging flag
  GetMacro(DebugFlag, bool);

  /// Print debugging messages if debugging is enabled
  void Debug(char *);
};


#endif


