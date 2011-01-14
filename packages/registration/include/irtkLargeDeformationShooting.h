/*=========================================================================
  Date      : $Date: 09.02.2010$
=========================================================================*/

#ifndef _IRTKLARGEDEFORMATIONSHOOTING_H
#define _IRTKLARGEDEFORMATIONSHOOTING_H

#include <irtkLargeDeformationSciCalcPack.h>

class EulerianShooting
{
  private:
  VectorField InvDiffeo, Diffeo, NablaI, NablaAdM, VelocityField,TempInvDiffeoLocal;
  VectorField AdjointVectorField,TempDiffeo,TempInvDiffeo,TempDiffeoLocal,TempVectorField;
  ScalarField ImTemplate, ImTarget, InitialMomentum, Image, Momentum, AdjointMomentum, AdjointImage,InputInitialMomentum;
  ScalarField GradientMomentum, TempAdMomentum, TempAdImage, OptimizedMomentum,TempScalarField,TempScalarField3;
  VectorField *TempInvDiffeoLocalTable;
  VectorField *TempDiffeoLocalTable;
  int NX;
  int NY;
  int NZ;
  int NT;
  float DeltaX;
  // number of iterations in the shooting method
  int IterationNumber;
  // time step of one iteration
  float DeltaTimeSubdiv;
  // Generic name for the limiter that can be SuperBee, UpWind or MinMod
  float (*Limiter)(float,float); 
  // Cost to optimize
  float Cost;
  float Energy;
  float MaxVectorField;
  public:

  // Constructor and Destructor
  EulerianShooting();
  ~EulerianShooting();


  void AllocateVariablesShooting();
  void InitializeAdjointVariables();
  void InitializeVariables();
  void ReadAndTreatInputImages();
  void ComputeVelocityField();
  void ComputeVelocityField(ScalarField *Momentum,VectorField * NablaI);
  void ComputeAdjointVectorField();
  void SchemeStep();
  void SchemeStep(VectorField * TempInvDiffeoLoc, VectorField * TempDiffeoLoc,VectorField * Output1, VectorField * Output2, int t1=0, int t2=0);
  void RungeKutta();
  void Scheme();
  void Shooting();
  void ShootingShow();
  void GradientDescent(int, float);
  void Gradient();
  void SaveResult();

  float SimilarityMeasure();
  static inline float MinModLimiter(float a, float b)
  {
    if (a*b>0.0)
    {
      if (a<b){return a;}
      else {return b;}
    }
    return 0.0;
  }
  static inline float UpWindLimiter(float a,float b){return 0.0;}

  static float SuperBeeLimiter(float , float);
  void Run();
  /// private parameters
  
  //fft convolver (to smooth the images)
  FFTconvolver3D FFTconvolver;
  
  /// public Parameters
  char SourceImageName[256];
  char TargetImageName[256];
  char InputInitialMomentumName[256];
  // Number of times from time 0 to time T
  int NbTimes;
  int NbIter;

  // weight of the norm in the cost function
  float alpha;
  
  // Number that represents the choice of the scheme in time (1 for RungeKutta third order otherwise simple Euler)
  int indicatorRungeKutta;
  double weight1,sigmaX1,sigmaY1,sigmaZ1;
  double weight2,sigmaX2,sigmaY2,sigmaZ2;
  double weight3,sigmaX3,sigmaY3,sigmaZ3;
  double weight4,sigmaX4,sigmaY4,sigmaZ4;
  int GreyLevAlign, Margin;
  float GLA_Padding_Src;         //if grey level alignment: padding value for the source image
  float GLA_Padding_Trg;         //if grey level alignment: padding value for the target image
  int indicatorInitialMomentum;
  int OutIniMoTxt;
  int OutVeloField;
  int OutDistEnSim;
  int OutDeformation;
  // size of the update in the GradientDescent
  float MaxUpdate;
  // Number that stands for the choice of the spatial scheme: 
  // 0 for UpWind, 2 for SuperBee, default is MinMod
  int indicatorLimiter;
  char PrefixOutputs[256];       //Prefix of the files containing the final velocity field
};

#endif
