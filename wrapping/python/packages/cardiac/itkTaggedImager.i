extern class itkTaggedImager
{
 public:
  /** The different tag pattern directions. */
  enum TagPatternDirection {TAG_DIRECTION_XY, TAG_DIRECTION_XZ,
                            TAG_DIRECTION_YZ, TAG_DIRECTION_XYZ};

  /** Constructor. Voxel sizes are in cm. */
  itkTaggedImager(double dEchoTime = 0.03, double dPulseRepetitionTime = 10,
    double dkx = 8, double dky = 8, double dkz = 0, double dTipAngle = 45,
    int nImageSizeX = 256, int nImageSizeY = 256, int nImageSizeZ = 10,
    double dVoxelSizeX = 0.05, double dVoxelSizeY = 0.05,
    double dVoxelSizeZ = 0.1, int nPaddingValue = -1,
    TagPatternDirection enumDirection = TAG_DIRECTION_XY);

  /** Destructor. */
  virtual ~itkTaggedImager();

  /** Intializes the imager. */
  void Initialize(double dEchoTime, double dPulseRepetitionTime, double dkx,
    double dky, double dkz, double dTipAngle, int nImageSizeX, int nImageSizeY,
    int nImageSizeZ, double dVoxelSizeX, double dVoxelSizeY,
    double dVoxelSizeZ, int nPaddingValue, TagPatternDirection enumDirection);

  /** Images an LV model at time nTimeInstant. */
  itkGreyImage Image(itkLVModel& lvModel, int nTimeInstant,
    bool outputDisplacements, ofstream* pOutputFile);

  /** Stores the displacement at each point in the myocardium in a set of
      real images. */
  void OutputDisplacements(itkLVModel& lvModel, int nTimeInstant,
    itkGreyImage& firstImage, itkRealImage& xDisps, itkRealImage& yDisps, itkRealImage& zDisps);

  /** Calculates the Jacobian of the transformation used by the LV model. */
  void CalculateJacobian(itkLVModel& lvModel, int nTimeInstant, int x, int y,
    int z, itkMatrix& jacobian);

  /** Stores the Jacobian at each point in the myocardium in a set of real
      images. */
  void OutputJacobians(itkLVModel& lvModel, int nTimeInstant, itkGreyImage& firstImage, itkRealImage& j00, itkRealImage& j01, itkRealImage& j02, itkRealImage& j10, itkRealImage& j11, itkRealImage& j12, itkRealImage& j20, itkRealImage& j21, itkRealImage& j22);

  /** Calculates the principal strain. */
  void CalculatePrincipalStrain(itkLVModel& lvModel, int nTimeInstant, int x,
    int y, int z, itkMatrix& principalStrainDirections,
    itkVector& principalStrains);

  /** Calculates the diagonal components of the strain tensor in Cartesian
      coordinates. */
  void CalculateStrainInCartesianCoordinates(itkLVModel& lvModel,
    int nTimeInstant, int x, int y, int z, double& ex, double& ey, double& ez);

  /** Calculates the components of the 3x3 strain tensor in Cartesian
      coordinates. */
  void CalculateStrainInCartesianCoordinates(itkLVModel& lvModel,
    int nTimeInstant, int x, int y, int z, itkMatrix& strainTensor);

  /** Calculates the specfied component of the strain tensor and returns the
      values in strainImage. */
  void CalculateStrainInCylindricalCoordinates(itkLVModel& lvModel,
    int nTimeInstant, itkSegmentedLVImage& segmentedLVImage,
    itkSegmentedLVImage::CylindricalStrainComponent component,
    itkRealImage& strainImage);

  /** Gets the echo time. */
  double GetEchoTime() const;

  /** Updates the scan image. */
  void UpdateScanImage();

  /** Gets the pulse repetition time. */
  double GetPulseRepetitionTime() const;

  /** Gets the x, y, z dimensions of the imaging volume. */
  void GetImageSize(int& nX, int& nY, int& nZ) const;

  /** Gets the voxel size. Returns value in cm. */
  void GetVoxelSize(double& dX, double& dY, double& dZ) const;

  /** Converts from image coordinates to world coordinates. */
  void ImageToWorld(double nX, double nY, double nZ, double& dX, double& dY,
    double& dZ) const;

  /** Converts from world coordinates to image coordiantes. */
  void WorldToImage(double dX, double dY, double dZ, double& nX, double& nY,
    double& nZ) const;

  /** Gets the spatial frequencies. */
  void GetSpatialFrequencies(double& dkx, double& dky, double& dkz) const;

  /** Gets the tip angle. */
  double GetTipAngle() const;
  
  /** Sets the direction of the tag pattern. */
  void SetTagPatternDirection(TagPatternDirection enumDirection);

  /** Gets the direction of the tag pattern. */
  TagPatternDirection GetTagPatternDirection(void) const;

  /** Sets the imaging volume origin. */
  void SetImagingVolumeOrigin(double x, double y, double z);

  /** Returns the imaging volume origin. */
  void GetImagingVolumeOrigin(double& x, double& y, double& z);

  /** Sets the imaging volume orientation. */
  void SetImagingVolumeOrientation(double xx, double xy, double xz, double yx,
    double yy, double yz);

  /** Gets the imaging volume orientation. */
  void GetImagingVolumeOrientation(double& xx, double& xy, double& xz,
    double& yx, double& yy, double& yz);
  
  /** Gets the padding value. */
  int GetPaddingValue() const;

  /** Sets parameters for the Gaussian noise which will be added to the image.
      Points not in the myocardium are not affected by the noise.
      If dSigma = 0 no noise is added. */
  void SetNoiseParameters(double dMean, double dSigma);

  /** Returns the noise parameters. */
  void GetNoiseParameters(double& dMean, double& dSigma) const;
};
