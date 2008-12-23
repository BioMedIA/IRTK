template <class VoxelType> extern class itkGaussianBlurring
{
public:
  itkGaussianBlurring(double);
  ~itkGaussianBlurring();
  virtual void SetInput(itkGenericImage<VoxelType>* );
  virtual void SetOutput(itkGenericImage<VoxelType>* );
  virtual void Run();
  void SetSigma(double);
  double GetSigma();
};
