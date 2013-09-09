/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/


/*
 * 	This class implements the Paper:
 * 	Neonatal brain segmentation using second order neighboorhood information, Ledig C. et al
 * 	PaPI Workshop 2012 in conjunction with MICCAI 2012, Nice *
 *
 */

#ifndef IRTKEMCLASSIFICATION2NDORDERMRF_H_
#define IRTKEMCLASSIFICATION2NDORDERMRF_H_

#include <irtkEMClassification.h>
#include <irtkPolynomialBiasField.h>
#include <irtkImage.h>
#include <irtkProbabilisticAtlas.h>
#include <irtkBiasField.h>
#include <irtkBiasCorrectionMask.h>
#include <map>
#include <irtkProbabilisticAtlas.h>

class irtkEMClassification2ndOrderMRF : public irtkEMClassification
{
protected:
  double _mrf_weight_distant;
  double _mrf_weight_adjacent;
  double _relax_factor;
  double _relax_factor_2nd;

  bool _isLogTransformed;

  /// local weights for MRF field
  irtkRealImage _MRF_weights;

  /// Uncorrected image
  irtkRealImage _uncorrected;

  /// Bias field correction filter
  irtkBiasCorrectionMask _biascorrection;

  /// Bias field
  irtkBiasField *_biasfield;

  // connectivity matrix (undirected->symmetric) (0 same class, 1 close class, 2 distant class)
  irtkMatrix _connectivity;

  vector< vector< pair< pair<int,int>, double > > > _connectivity_2nd_order;
  bool _has_MRF_2nd_order;

  bool _has_background;

  map<int,int> pv_classes;
  vector< pair<int, int> > pv_connections;
  vector<double> pv_fc;

private:
  bool isPVclass(int pvclass);
  double getMRFenergy(int index, int tissue);
  double getMRFenergy_2nd_order(int index, int tissue);

  double getTau(int index, int tissue);

public:
  /// Constructor
  irtkEMClassification2ndOrderMRF();

  /**
      * Constructor
      * @param noTissues number of different tissue types to be segmented
      * @param atlas pointer to atlases, number of given atlases should conincide with noTissues
      * @param background pointer to background image
  */
  irtkEMClassification2ndOrderMRF(int noTissues, irtkRealImage **atlas, irtkRealImage *background);
  irtkEMClassification2ndOrderMRF(int noTissues, irtkRealImage **atlas);

  void SetLogTransformed(bool);
  double LogLikelihood();
  void EStep();

  void GetBiasField(irtkRealImage &image);

  void SetMRFWeights(double adjacent, double distant);

  void SetRelaxationFactor(double relax_factor);

  void SetRelaxationFactor_2nd(double relax_factor);

  void SetMRF_2nd_order(vector< vector< pair< pair<int,int>, double > > > &);

  /// sets manually the distribution parameters
  void SetParameters( double mu, double var, int number);

  /// removes the PV classes!
  void removePVclasses(double threshold = 0.1);

  /// construct segmentation with existing PV classes!
  void ConstructSegmentationPV(irtkRealImage &segmentation);

  void ConstructSegmentation(irtkRealImage &segmentation);

  void RefineSegmentation(int);

  /// add partial volume between classes classA and classB
  int AddPartialVolumeClass(int classA, int classB);

  /// estimate probabilities using 1st order MRF only
  void EStepMRF(void);

  /// estimate probabilities using 1st and 2nd order MRF
  void EStepMRF_2nd_order(void);

  /// relax priors
  void RStep(void);

  /// relax priors based on 2nd order information
  void RStep_2nd_order(void);

  /// Estimates bias field
  virtual void BStep();

  // corrected M step for partial volume classes
  virtual void MStepPV();

  /// Set image
  virtual void SetInput(const irtkRealImage &, const irtkMatrix &);

  /// Set bias field
  virtual void SetBiasField(irtkBiasField *);

  /// Compute the bias corrected image
  virtual void GetBiasCorrectedImage(irtkRealImage &);


};

#endif /* irtkEMClassification2ndOrderMRF_H_ */
