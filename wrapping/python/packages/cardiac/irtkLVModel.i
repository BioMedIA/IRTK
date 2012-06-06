extern class itkLVModel
{
public:
  /** 13 k-parameters are needed to describe the model. */
  static const int NO_K_PARAMS = 13;

  /** MINETA and MAXETA describe the angular extent of the LV model. */ 
  static const double MINETA;
  static const double MAXETA;

  /** Constructor. */
  itkLVModel(int nTimeInstants = 20, double dLengthTimeInterval = 1,
    double dSpinDensity = 300.0, double dLongRelaxationTime = 0.6,
    double dTransverseRelaxationTime = 0.10, double dDelta = 4.00,
    double dLambdai = 0.35, double dLambdao = 0.55,
    const double* ppdKParams[NO_K_PARAMS] = NULL);

  /** Initializes the model. */
  void Initialize(int nTimeInstants, double dLengthTimeInterval,
    double dSpinDensity, double dLongRelaxationTime,
    double dTransverseRelaxationTime, double dDelta, double dLambdai,
    double dLambdao, const double* ppdKParams[NO_K_PARAMS]);

  /** Gets the number of time instants. */
  int GetNoTimeInstants() const;

  /** Gets the duration of the modelling. */
  double GetLengthTimeInterval() const;

  /** Gets the spin density of the model. */
  double GetSpinDensity() const;

  /** Gets the longitudianl relaxation time. */
  double GetLongitudinalRelaxationTime() const;

  /** Gets the transverse relaxation time. */
  double GetTransverseRelaxationTime() const;

  /** Returns true if a point is in the model at time t = 0. */
  bool IsPointInModel(const itkMatrix& matPosition) const;

  /** Gets the model dimensions. */
  void GetModelDimensions(double& dDelta, double& dLambdai, double& dLambdao);

  /** Converts from prolate spheroidal coordinates to Cartesian coordinates. */
  void ProlateSpheroidalToCartesian(double lambda, double eta, double phi,
                                    double& OUTPUT, double& OUTPUT, double& OUTPUT);

  /** Converts from Cartesian coordinates to prolate spheroidal coordinates. */
  void CartesianToProlateSpheroidal(double x, double y, double z,
                                    double& OUTPUT, double& OUTPUT, double& OUTPUT);

  /** Calculates the transformed position of a point in the LV model. */
  itkMatrix CalculateTransformationPosition(const itkMatrix& refPosition,
    int nTimeInstant);

  /** Calculates the displacement of a point in the LV model at time
      nTimeInstant from its position at time t = 0. */
  itkMatrix CalculateDisplacement(const itkMatrix& position, int nTimeInstant);

  /** Calculates the reference map, which transforms points in the LV model at
      time t to their positions at time 0. */
  itkMatrix CalculateReferencePosition(const itkMatrix& matPosition,
    int nTimeInstant);

  /** Sets a k-parameter value.
      \param k The k-parameter to set.
      \param time The time frame.
      \param val The new value. */
  void SetKParameterValue(int k, int time, double val);

  /** Returns a k-parameter value.
      \param k The k-parameter to get.
      \param time The time frame. */
  double GetKParameterValue(int k, int time);

  /** Initializes the model by computing the position independent matrices.
      This should be called after any of the k-parameters have been modified. */
  void InitializeMatrices();
};
