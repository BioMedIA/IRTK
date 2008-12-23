extern class itkMultiFrameLatticeFreeFormTransformation : public itkTransformation
{
public:
  itkMultiFrameLatticeFreeFormTransformation();
  virtual ~itkLatticeFreeFormTransformation();
  virtual int NumberOfDOFs() const;
  virtual double Get(int p) const;
  virtual void Put(int p, double val);
  virtual bool IsTransformDefinedAtPoint(double x, double y, double z);
  virtual void Transform(double& INOUT, double& INOUT, double& INOUT);
  virtual void Jacobian(double x, double y, double z, itkMatrix& jac);
  virtual void LocalJacobian(double x, double y, double z, itkMatrix& jac);
  virtual void GlobalJacobian(double x, double y, double z, itkMatrix& jac);
  virtual Bool IsIdentity();
  virtual void Print();
  virtual char* NameOfClass();
  void UpdateVertexPositions();
  unsigned int GetMaximumLevel() const;
  virtual void SetTimeFrame(unsigned int time);
  virtual unsigned int GetTimeFrame() const;
  virtual unsigned int GetNumberOfTimeFrames() const;
  void InitializeCurrentTimeFrame();
%extend
{
  void Read(char* pFileName)
  {
    self->itkTransformation::Read(pFileName);  
  }
}
};
