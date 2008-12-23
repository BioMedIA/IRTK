extern class itkQuaternionTransformation : public itkTransformation
{
public:
  %apply double DOUBLE3INPUT[3] {double originInput[3], double alphaAxisInput[3], double betaAxisInput[3]};
  %apply double DOUBLE3OUTPUT[3] {double originOutput[3], double alphaAxisOutput[3], double betaAxisOutput[3], double gammaAxisOutput[3]};
  
  itkQuaternionTransformation();
  void SetAlpha(double angle);
  void SetBeta(double angle);
  void SetGamma(double angle);
  double GetAlpha() const;
  double GetBeta() const;
  double GetGamma() const;
  void SetAxes(double alphaAxisInput[3], double betaAxisInput[3]);
  void GetAxes(double alphaAxisOutput[3], double betaAxisOutput[3], double gammaAxisOutput[3]) const;
  void GetRotatedAxes(double alphaAxisOutput[3], double betaAxisOutput[3], double gammaAxisOutput[3]) const;
  void SetOrigin(double originInput[3]);
  void GetOrigin(double originOutput[3]) const;
  const itkMatrix& GetMatrix() const;
  void UpdateMatrix();
  virtual int NumberOfDOFs() const;
  virtual double Get(int param) const;
  virtual void Put(int param, double val);
  virtual void Transform(double& INOUT, double& INOUT, double& INOUT);
  virtual void Invert();
  virtual void Jacobian(double x, double y, double z, itkMatrix& jac);
  virtual double Jacobian(double x, double y, double z);
  virtual void LocalJacobian(double x, double y, double z, itkMatrix& jac);
  virtual void GlobalJacobian(double x, double y, double z, itkMatrix& jac);
  static int CheckHeader(char* pFileName);
  virtual Bool IsIdentity();
  virtual void Print();
  virtual char* NameOfClass();
};
