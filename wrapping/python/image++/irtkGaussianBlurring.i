template <class VoxelType> extern class irtkGaussianBlurring
{
public:
  irtkGaussianBlurring(double);
  ~irtkGaussianBlurring();
  virtual void SetInput(irtkGenericImage<VoxelType>* );
  virtual void SetOutput(irtkGenericImage<VoxelType>* );
  virtual void Run();
  void SetSigma(double);
  double GetSigma();
};
