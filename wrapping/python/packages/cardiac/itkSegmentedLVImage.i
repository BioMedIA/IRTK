extern class itkSegmentedLVImage
{
 public:
  /** Constructor. */
  itkSegmentedLVImage();

  /** Destructor. */
  virtual ~itkSegmentedLVImage();

  /** Sets the segmented LV image. */
  void SetSegmentedImage(itkGreyImage* pSegmentedImage, itkGreyPixel nPaddingValue);

  /** Gets the segmented LV image. */
  itkGreyImage* GetSegmentedImage(void) const;
  
  /** Gets the padding value. */
  itkGreyPixel GetPaddingValue(void) const;

  /** Returns true if a particular slice in the segmented LV image contains
   *  the myocardium.
   */
  bool IsMyocardiumInSlice(int z) const;

  /** Calculates the axis transformation for the segmented image. This
   *  transformation moves from a coordinate system, where the z-axis is
   *  aligned with a line passing from the center of the base to the center of
   *  the apex, to the world coordinate system. The transformation returned
   *  must be deleted by the user.
   */
  itkRigidTransformation* CalculateRigidAxisTransformation(void) const;
  
  /** Same as above function but a matrix is returned. */
  void CalculateRigidAxisTransformation(itkMatrix& matrix) const;

  /** Computes the rigid axis transformation for the segmented image. This
      transformation moves from a coordinate system, where the z-axis is
      aligned with a line passing from the center of the base to the center of
      the apex, to the world coordinate system. The transformation returned
      must be deleted by the user.
      \param axis The direction of the vector joining the center of the LV to
      the lateral part of the LV.
  */
  itkRigidTransformation* CalculateRigidAxisTransformation(itkPoint& axis) const;
  /** Same as above function but a matrix is returned. */
  void CalculateRigidAxisTransformation(itkMatrix& totalMatrix, itkPoint& xAxis) const;

  /** Calculates the rigid axis transformation needed for use with an
      itkLVTransformation4D transformation. */
  void CalculateLV4DRigidAxisTransformation(itkMatrix& matrix) const;
  
  /** Calculates the axis transformation for the segmented image. This
   *  transformation moves from a coordinate system, where the z-axis is
   *  aligned with a line passing from the center of the base to the center of
   *  the apex, to the world coordinate system. The transformation returned
   *  must be deleted by the user.
   */
  itkAffineTransformation* CalculateAffineAxisTransformation(void) const;

  enum CartesianStrainComponent {
    CAR_XX, CAR_XY, CAR_XZ,
    CAR_YX, CAR_YY, CAR_YZ,
    CAR_ZX, CAR_ZY, CAR_ZZ
  };
  /** Calculates the specified component of the strain tensor for all
      myocardial points in the segmented image of the LV. The strain
      is returned in strain image. */
  void CalculateStrainInCartesianCoordinateSystem(itkTransformation& transformation, CartesianStrainComponent component, itkRealImage& strainImage);
  
  /** Calculates strain tensor in Cartesian coordinate system. */
  void CalculateStrainInCartesianCoordinateSystem(itkTransformation& transformation, int x, int y, int z, itkMatrix& strain);

  /** Computes the Lagrangian strain tensor at a point. */
  void CalculateLagrangianStrain(itkTransformation& transformation, double x, double y, double z, itkMatrix& strain);

  /** Computes the Lagrangian strain tensor in the cylindrical coordinate
      system.
      \param transformation Is the transformation.
      \param wx The x-position (in world coordinates).
      \param wy The y-position (in world coordinates).
      \param wz The z-position (in world coordinates).
      \param pAxisTransformation The axis transformation.
      \param pInverseAxisTransformation The inverse axis transformation.
      \param strain The strain. */
  void CalculateLagrangianStrainInCylindricalCoordinateSystem(itkTransformation& transformation, double wx, double wy, double wz, itkRigidTransformation* pAxisTransformation, itkRigidTransformation* pInverseAxisTransformation, itkMatrix& strain);

  /** Computes the radial, circumferential, and longitudinal directions in
      the world coordinate system.
      \param wx The x-coordinate of the point.
      \param wy The y-coordinate of the point.
      \param wz The z-coordinate of the point.
      \param pAxisTransformation The axis transformation.
      \param pInverseAxisTransformation The inverse axis transformation.
      \param rad The radial direction.
      \param circ The circumferential direction.
      \param longi The longitudinal direction. */
  void CalculateRCLDirections(double wx, double wy, double wz, itkRigidTransformation* pAxisTransformation, itkRigidTransformation* pInverseAxisTransformation, itkPoint& rad, itkPoint& circ, itkPoint& longi);
  
  /** Calculates strain tensor in Cartesian coordinate system. */
  void CalculateStrainInCartesianCoordinateSystem(itkMultiLevelCylindricalFreeFormTransformation& transformation, int x, int y, int z, itkMatrix& strain);

  enum CylindricalStrainComponent {
    CYL_RR, CYL_RT, CYL_RZ,
    CYL_TR, CYL_TT, CYL_TZ,
    CYL_ZR, CYL_ZT, CYL_ZZ
  };
  
  /** Calculates the specified component of the strain tensor for all
      myocardial points in the segmented image of the LV. The strain
      is returned in strain image. */
  void CalculateStrainInCylindricalCoordinateSystem(itkTransformation& transformation, CylindricalStrainComponent component, itkRealImage& strainImage);
  
  /** Calculates the specified component of the strain tensor for all
      myocardial points in the segmented image of the LV for cylindrical
      free-form transformations. The strain is returned in the strain image.
  */
  void CalculateStrainInCylindricalCoordinateSystem(itkMultiLevelCylindricalFreeFormTransformation& transformation, CylindricalStrainComponent component, itkRealImage& strainImage);
  
  /** Calculates the principal strain. */
  void CalculatePrincipalStrain(double wx, double wy, double wz, itkTransformation& transformation, itkMatrix& principalStrainDirections, itkVector& principalStrains);
};
