/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkShapeBasedInterpolateImageFunction.cc 8 2009-03-02 16:12:58Z dr $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-03-02 16:12:58 +0000 (一, 02 三月 2009) $
  Version   : $Revision: 8 $
  Changes   : $Author: dr $

  =========================================================================*/

#include <irtkImage.h>

#include <irtkImageFunction.h>

#include <irtkResampling.h>

#ifdef HAS_CONTRIB

#include <irtkEuclideanDistanceTransform.h>

#endif

irtkShapeBasedInterpolateImageFunction::irtkShapeBasedInterpolateImageFunction()
{
}

irtkShapeBasedInterpolateImageFunction::~irtkShapeBasedInterpolateImageFunction(void)
{
}

const char *irtkShapeBasedInterpolateImageFunction::NameOfClass()
{
  return "irtkShapeBasedInterpolateImageFunction";
}

#ifdef HAS_CONTRIB

void irtkShapeBasedInterpolateImageFunction::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

  double xsize, ysize, zsize, size;
  int new_x, new_y, new_z, x, y, z, t,labelcount;
  double xaxis[3], yaxis[3], zaxis[3];
  double new_xsize, new_ysize, new_zsize;
  double old_xsize, old_ysize, old_zsize;
  double min,max,current,sum,sumcount;

  _dmap = irtkRealImage(_input->GetImageAttributes());
  _tinput = irtkRealImage(_input->GetImageAttributes());

  // Initialize _rinput _rdmap
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  size = xsize;
  size = (size < ysize) ? size : ysize;
  size = (size < zsize) ? size : zsize;
  if(size > 1) size = 1;
  cerr << "Create Images with isotropic voxel size (in mm): "<< size << endl;

  // Create _rinput _rdmap

  // Determine the old dimensions of the image
  _input->GetPixelSize(&old_xsize, &old_ysize, &old_zsize);

  // Determine the new dimensions of the image
  new_x = int(_input->GetX() * old_xsize / size);
  new_y = int(_input->GetY() * old_ysize / size);
  new_z = int(_input->GetZ() * old_zsize / size);

  // Determine the new voxel dimensions
  if (new_x < 1) {
    new_x     =  1;
    new_xsize =  old_xsize;
  } else {
    new_xsize = size;
  }
  if (new_y < 1) {
    new_y     =  1;
    new_ysize =  old_ysize;
  } else {
    new_ysize = size;
  }
  if (new_z < 1) {
    new_z     =  1;
    new_zsize =  old_zsize;
  } else {
    new_zsize = size;
  }

  // Allocate new image
  _rinput = irtkRealImage(new_x, new_y, new_z, _input->GetT());
  _rdmap = irtkRealImage(new_x, new_y, new_z, _input->GetT());
  _rcdmap = irtkRealImage(new_x, new_y, new_z, _input->GetT());

  // Set new voxel size
  _rinput.PutPixelSize(new_xsize, new_ysize, new_zsize);
  _rdmap.PutPixelSize(new_xsize, new_ysize, new_zsize);
  _rcdmap.PutPixelSize(new_xsize, new_ysize, new_zsize);

  // Set new orientation
  _input ->GetOrientation(xaxis, yaxis, zaxis);
  _rinput.PutOrientation(xaxis, yaxis, zaxis);
  _rdmap.PutOrientation(xaxis, yaxis, zaxis);
  _rcdmap.PutOrientation(xaxis, yaxis, zaxis);

  // Set new origin
  _rinput.PutOrigin(_input->GetOrigin());
  _rdmap.PutOrigin(_input->GetOrigin());
  _rcdmap.PutOrigin(_input->GetOrigin());

  // For every intensity value
  _input->GetMinMaxAsDouble(&min,&max);

  labelcount = 0;
  for(current = min; current <= max; current ++){
    // Threshold
    sumcount = 0;
    for (t = 0; t < _tinput.GetT(); t++){
      for (z = 0; z < _tinput.GetZ(); z++) {
	for (y = 0; y < _tinput.GetY(); y++) {
	  for (x = 0; x < _tinput.GetX(); x++) {
	    if (_input->GetAsDouble(x, y, z, t) < current
		//|| _input->GetAsDouble(x, y, z, t) > current
		) {
	      _tinput(x, y, z, t) = 0;
	    } else {
	      _tinput(x, y, z, t) = 1;
	    }
	    if(_input->GetAsDouble(x, y, z, t) >= current
	       && _input->GetAsDouble(x, y, z, t) < current + 1){
	      sumcount ++;
	    }
	  }
	}
      }
    }
	  
    if(round(current)%20 == 0){
      if(max < current + 19){
	cout << "Doing outside DT for value >= : "<< current << " to " << max << endl;
	cout << "Doing inside DT for value >= : "<< current << " to " << max << endl;
      }else{
	cout << "Doing outside DT for value >= : "<< current << " to " << current+19 << endl;
	cout << "Doing inside DT for value >= : "<< current << " to " << current+19 << endl;
      }
    }

    if(sumcount > 0){
      labelcount ++;
      // Dmap _tinput to _dmap
      {
	irtkRealImage inputA, inputB, outputA, outputB;

	// Default mode
	irtkEuclideanDistanceTransform<irtkRealPixel> 
	  *edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
	  (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);

	// Threshold image
	inputA = _tinput;
	inputB = _tinput;
	for (t = 0; t < _tinput.GetT(); t++){
	  for (z = 0; z < _tinput.GetZ(); z++) {
	    for (y = 0; y < _tinput.GetY(); y++) {
	      for (x = 0; x < _tinput.GetX(); x++) {
		if (_tinput(x, y, z, t) > 0.5) {
		  inputA(x, y, z, t) = 1;
		  inputB(x, y, z, t) = 0;
		} else {
		  inputA(x, y, z, t) = 0;
		  inputB(x, y, z, t) = 1;
		}
	      }
	    }
	  }
	}

	// Calculate EDT
	edt->SetInput (& inputA);
	edt->SetOutput(&outputA);
	edt->Run();		  
	edt->SetInput (& inputB);
	edt->SetOutput(&outputB);
	edt->Run();
	for (t = 0 ; t < _tinput.GetT(); t++){
	  for (z = 0; z < _tinput.GetZ(); z++) {
	    for (y = 0; y < _tinput.GetY(); y++) {
	      for (x = 0; x < _tinput.GetX(); x++) {
		_dmap(x, y, z, t)  = sqrt(outputA(x, y, z, t)) - sqrt(outputB(x, y, z, t));
	      }
	    }
	  }
	}
	//fix the result to better visiualization and cspline interpolation
	/*_dmap.GetMinMaxAsDouble(&dmin,&dmax);
	  if( abs(dmin) > abs(dmax))
	  dmin = abs(dmax);
	  for (t = 0; t < _tinput.GetT(); t++){
	  for (z = 0; z < _tinput.GetZ(); z++) {
	  for (y = 0; y < _tinput.GetY(); y++){
	  for (x = 0; x < _tinput.GetX(); x++){
	  if(_dmap.GetAsDouble(x,y,z,t) > abs(dmin)){
	  _dmap(x, y, z, t) = abs(dmin);
	  }
	  }
	  }			
	  }
	  }*/
	//_dmap.Write("testoutput0.nii");
      }

      // Linear Interpolate Dmap _dmap to _rdmap
      {
	double i,j,k,l;
	irtkLinearInterpolateImageFunction interpolator;
	interpolator.SetInput(&_dmap);
	interpolator.Initialize();

	for (t = 0 ; t < _rdmap.GetT(); t++){
	  for (z = 0; z < _rdmap.GetZ(); z++) {
	    for (y = 0; y < _rdmap.GetY(); y++) {
	      for (x = 0; x < _rdmap.GetX(); x++) {
		i = x; j = y; k = z; l = t;
		_rdmap.ImageToWorld(i,j,k);
		_dmap.WorldToImage(i,j,k);
		_rdmap.PutAsDouble(x,y,z,t,interpolator.Evaluate(i,j,k,l));
	      }
	    }
	  }
	}
	//_rdmap.Write("testoutput1.nii");
      }

      // Put value back to Resampled Image _rinput < 0 if _rinput == 0 let the neighbor vote.
      {
	for (t = 0 ; t < _rdmap.GetT(); t++){
	  for (z = 0; z < _rdmap.GetZ(); z++) {
	    for (y = 0; y < _rdmap.GetY(); y++) {
	      for (x = 0; x < _rdmap.GetX(); x++) {
		if(_rdmap.GetAsDouble(x,y,z,t) < 0){
		  _rinput.PutAsDouble(x,y,z,t,current);
		}else if(_rdmap.GetAsDouble(x,y,z,t) <= 0){
		  sum = 0; sumcount = 0;
		  if (x > 0) {
		    sum += _rdmap.GetAsDouble(x-1,y,z,t); 
		    sumcount ++;
		  }
		  if (x < _rinput.GetX() - 1) {
		    sum += _rdmap.GetAsDouble(x+1,y,z,t); 
		    sumcount ++;
		  }
		  if (y > 0) {
		    sum += _rdmap.GetAsDouble(x,y-1,z,t); 
		    sumcount ++;
		  } 
		  if (y < _rinput.GetY() - 1) {
		    sum += _rdmap.GetAsDouble(x,y+1,z,t); 
		    sumcount++;
		  }
		  if (z > 0) {
		    sum += _rdmap.GetAsDouble(x,y,z-1,t); 
		    sumcount++;
		  }
		  if (z < _rinput.GetZ() - 1) {
		    sum += _rdmap.GetAsDouble(x,y,z+1,t); 
		    sumcount++;
		  }
		  sum = sum/sumcount;
		  if(sum <= 0){
		    _rinput.PutAsDouble(x,y,z,t,current);
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // End for
  }
  if(labelcount > 3 && labelcount < 50){
    //_rinput.Write("beforerefine.nii");
    this->Refine();
  }
  // Compute image domain
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  this->_rx = this->_rinput.GetX();
  this->_ry = this->_rinput.GetY();
  this->_rz = this->_rinput.GetZ();

  // Compute domain on which the linear interpolation is defined
  this->_x1 = 0;
  this->_y1 = 0;
  this->_z1 = 0;
  this->_x2 = this->_rinput.GetX() - 1;
  this->_y2 = this->_rinput.GetY() - 1;
  this->_z2 = this->_rinput.GetZ() - 1;

  // Calculate offsets for fast pixel access
  this->_offset1 = 0;
  this->_offset2 = 1;
  this->_offset3 = this->_rinput.GetX();
  this->_offset4 = this->_rinput.GetX()+1;
  this->_offset5 = this->_rinput.GetX()*this->_rinput.GetY();
  this->_offset6 = this->_rinput.GetX()*this->_rinput.GetY()+1;
  this->_offset7 = this->_rinput.GetX()*this->_rinput.GetY()+this->_rinput.GetX();
  this->_offset8 = this->_rinput.GetX()*this->_rinput.GetY()+this->_rinput.GetX()+1;
}

void irtkShapeBasedInterpolateImageFunction::Refine()
{

  int x, y, z, t;
  double min,max,current,sum,sumcount,dcurrent;

  // Initialization
  sum = 0;

  // For every intensity value
  _input->GetMinMaxAsDouble(&min,&max);

  // Initialize _rcdmap
  for (t = 0 ; t < _rcdmap.GetT(); t++){
    for (z = 0; z < _rcdmap.GetZ(); z++) {
      for (y = 0; y < _rcdmap.GetY(); y++) {
	for (x = 0; x < _rcdmap.GetX(); x++) {
	  _rcdmap.PutAsDouble(x,y,z,t,30);
	}
      }
    }
  }

  for(current = min; current <= max; current ++){
    // Threshold
    sumcount = 0;
    for (t = 0; t < _tinput.GetT(); t++){
      for (z = 0; z < _tinput.GetZ(); z++) {
	for (y = 0; y < _tinput.GetY(); y++) {
	  for (x = 0; x < _tinput.GetX(); x++) {
	    if ((_input->GetAsDouble(x, y, z, t) < current)
		|| (_input->GetAsDouble(x, y, z, t) >= (current + 1))
		) {
	      _tinput(x, y, z, t) = 0;
	    } else {
	      _tinput(x, y, z, t) = 1;
	      sumcount ++;
	    }
	  }
	}
      }
    }

    // Calculate EDT
    if(round(current)%20 == 0){
      if(max < current + 19){
	cout << "Doing outside DT for value == : "<< current << " to " << max << endl;
	cout << "Doing inside DT for value == : "<< current << " to " << max << endl;
      }else{
	cout << "Doing outside DT for value == : "<< current << " to " << current+19 << endl;
	cout << "Doing inside DT for value == : "<< current << " to " << current+19 << endl;
      }
    }

    if(sumcount > 0){
      // Dmap _tinput to _dmap
      {
	irtkRealImage inputA, inputB, outputA, outputB;

	// Default mode
	irtkEuclideanDistanceTransform<irtkRealPixel> 
	  *edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
	  (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);
		  
	// Threshold image
	inputA = _tinput;
	inputB = _tinput;
	for (t = 0; t < _tinput.GetT(); t++){
	  for (z = 0; z < _tinput.GetZ(); z++) {
	    for (y = 0; y < _tinput.GetY(); y++) {
	      for (x = 0; x < _tinput.GetX(); x++) {
		if (_tinput(x, y, z, t) > 0.5) {
		  inputA(x, y, z, t) = 1;
		  inputB(x, y, z, t) = 0;
		} else {
		  inputA(x, y, z, t) = 0;
		  inputB(x, y, z, t) = 1;
		}
	      }
	    }
	  }
	}

	edt->SetInput (& inputA);
	edt->SetOutput(&outputA);
	edt->Run();		  
	edt->SetInput (& inputB);
	edt->SetOutput(&outputB);
	edt->Run();
	for (t = 0 ; t < _tinput.GetT(); t++){
	  for (z = 0; z < _tinput.GetZ(); z++) {
	    for (y = 0; y < _tinput.GetY(); y++) {
	      for (x = 0; x < _tinput.GetX(); x++) {
		_dmap(x, y, z, t)  = sqrt(outputA(x, y, z, t)) - sqrt(outputB(x, y, z, t));
	      }
	    }
	  }
	}
      }

      /*_dmap.GetMinMaxAsDouble(&dmin,&dmax);
	if( abs(dmin) > abs(dmax))
	dmin = abs(dmax);
	for (t = 0; t < _tinput.GetT(); t++){
	for (z = 0; z < _tinput.GetZ(); z++) {
	for (y = 0; y < _tinput.GetY(); y++){
	for (x = 0; x < _tinput.GetX(); x++){
	if(_dmap.GetAsDouble(x,y,z,t) > abs(dmin)){
	_dmap(x, y, z, t) = abs(dmin);
	}
	}
	}			
	}
	}*/

      // Linear Interpolate Dmap _dmap to _rdmap
      {
	double i,j,k,l;
	irtkLinearInterpolateImageFunction interpolator;
	interpolator.SetInput(&_dmap);
	interpolator.Initialize();

	for (t = 0 ; t < _rdmap.GetT(); t++){
	  for (z = 0; z < _rdmap.GetZ(); z++) {
	    for (y = 0; y < _rdmap.GetY(); y++) {
	      for (x = 0; x < _rdmap.GetX(); x++) {
		i = x; j = y; k = z; l = t;
		_rdmap.ImageToWorld(i,j,k);
		_dmap.WorldToImage(i,j,k);
		_rdmap.PutAsDouble(x,y,z,t,interpolator.Evaluate(i,j,k,l));
	      }
	    }
	  }
	}
	//_rdmap.Write("testoutput1.nii");
      }

      // Put value back to Resampled Image _rinput < 0 if _rinput == 0 let the neighbor vote.
      {
	for (t = 0 ; t < _rdmap.GetT(); t++){
	  for (z = 0; z < _rdmap.GetZ(); z++) {
	    for (y = 0; y < _rdmap.GetY(); y++) {
	      for (x = 0; x < _rdmap.GetX(); x++) {
		dcurrent = _rdmap.GetAsDouble(x,y,z,t);
		if(dcurrent < 0 
		   && current >= _rinput.GetAsDouble(x,y,z,t)
		   && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)){
		  //_rinput.PutAsDouble(x,y,z,t,current);
		  _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
		}else if(dcurrent <= 0 
			 && current >= _rinput.GetAsDouble(x,y,z,t)
			 && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)){
		  sum = 0; sumcount = 0;
		  if (x > 0) {
		    sum += _rdmap.GetAsDouble(x-1,y,z,t); 
		    sumcount ++;
		  }
		  if (x < _rinput.GetX() - 1) {
		    sum += _rdmap.GetAsDouble(x+1,y,z,t); 
		    sumcount ++;
		  }
		  if (y > 0) {
		    sum += _rdmap.GetAsDouble(x,y-1,z,t); 
		    sumcount ++;
		  } 
		  if (y < _rinput.GetY() - 1) {
		    sum += _rdmap.GetAsDouble(x,y+1,z,t); 
		    sumcount++;
		  }
		  if (z > 0) {
		    sum += _rdmap.GetAsDouble(x,y,z-1,t); 
		    sumcount++;
		  }
		  if (z < _rinput.GetZ() - 1) {
		    sum += _rdmap.GetAsDouble(x,y,z+1,t); 
		    sumcount++;
		  }
		  sum = sum/sumcount;
		  if(sum <= 0){
		    //_rinput.PutAsDouble(x,y,z,t,current);
		    _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
		  }
		}else if(dcurrent > 0 
			 && _rinput.GetAsDouble(x,y,z,t) <= current
			 && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)){
		  _rinput.PutAsDouble(x,y,z,t,current);
		  _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
		}else if(dcurrent >= 0 
			 && _rinput.GetAsDouble(x,y,z,t) <= current
			 && dcurrent < _rcdmap.GetAsDouble(x,y,z,t)){							  
		  if(sum > 0){
		    _rinput.PutAsDouble(x,y,z,t,current);
		    _rcdmap.PutAsDouble(x,y,z,t,dcurrent);
		  }
		}
	      }
	    }
	  }
	}
	//_rinput.Write("testoutput2.nii");
      }
    }

    // End for
  }
}

double irtkShapeBasedInterpolateImageFunction::EvaluateInside(double x, double y, double z, double t)
{
  int i, j, k, l;

  _input->ImageToWorld(x,y,z);
  _rinput.WorldToImage(x,y,z);

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  return this->_rinput.GetAsDouble(i, j, k, l);
}

double irtkShapeBasedInterpolateImageFunction::Evaluate(double x, double y, double z, double t)
{
  int i, j, k, l;

  _input->ImageToWorld(x,y,z);
  _rinput.WorldToImage(x,y,z);

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  if ((i < 0) || (i >= this->_rx) || (j < 0) || (j >= this->_ry) || (k < 0) || (k >= this->_rz)) {
    return this->_DefaultValue;
  } else {
    return this->_rinput.GetAsDouble(i, j, k, l);
  }
}

double irtkShapeBasedInterpolateImageFunction::EvaluateInsideLinear(double x, double y, double z, double t)
{
  int i, j, k;
  double t1, t2, u1, u2, v1, v2;

  _input->ImageToWorld(x,y,z);
  _rinput.WorldToImage(x,y,z);

  // Calculated integer coordinates
  i  = int(x);
  j  = int(y);
  k  = int(z);

  // Calculated fractional coordinates
  t1 = x - i;
  u1 = y - j;
  v1 = z - k;
  t2 = 1 - t1;
  u2 = 1 - u1;
  v2 = 1 - v1;

  // Get pointer to data
  switch (this->_input->GetScalarType()) {
  case IRTK_VOXEL_UNSIGNED_SHORT: {
    unsigned short *ptr = (unsigned short *)this->_rinput.GetScalarPointer(i, j, k, round(t));

    // Linear interpolation
    return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
		  u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
	    t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
		  u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

    break;
  }
  case IRTK_VOXEL_SHORT: {
    short *ptr = (short *)this->_rinput.GetScalarPointer(i, j, k, round(t));

    // Linear interpolation
    return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
		  u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
	    t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
		  u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

    break;
  }
  case IRTK_VOXEL_FLOAT: {
    float *ptr = (float *)this->_rinput.GetScalarPointer(i, j, k, round(t));

    // Linear interpolation
    return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
		  u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
	    t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
		  u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

    break;
  }
  case IRTK_VOXEL_DOUBLE: {
    double *ptr = (double *)this->_rinput.GetScalarPointer(i, j, k, round(t));

    // Linear interpolation
    return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
		  u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
	    t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
		  u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

    break;
  }
  default:
    cerr << "irtkLinearInterpolateImageFunction::EvaluateInside: Unknown scalar type" << endl;
    exit(1);
  }
}

double irtkShapeBasedInterpolateImageFunction::EvaluateLinear(double x, double y, double z, double t)
{
  int i, j, k, l, m, n;
  double val;

  _input->ImageToWorld(x,y,z);
  _rinput.WorldToImage(x,y,z);

  i = (int)floor(x);
  j = (int)floor(y);
  k = (int)floor(z);
  t = round(t);

  val = 0;
  for (l = i; l <= i+1; l++) {
    if ((l >= 0) && (l < this->_rx)) {
      for (m = j; m <= j+1; m++) {
        if ((m >= 0) && (m < this->_ry)) {
          for (n = k; n <= k+1; n++) {
            if ((n >= 0) && (n < this->_rz)) {
              val += (1 - fabs(l - x))*(1 - fabs(m - y))*(1 - fabs(n - z))*this->_rinput.GetAsDouble(l, m, n, t);
            }
          }
        }
      }
    }
  }
  switch (this->_input->GetScalarType()) {
  case IRTK_VOXEL_UNSIGNED_SHORT: {
    val = round(val);
    break;
  }
  case IRTK_VOXEL_SHORT: {
    val = round(val);
    break;
  }
  default: break;
  }
  return val;
}

#else

void irtkShapeBasedInterpolateImageFunction::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

  // Compute image domain
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  // Compute domain on which the linear interpolation is defined
  this->_x1 = -0.5;
  this->_y1 = -0.5;
  this->_z1 = -0.5;
  this->_x2 = this->_input->GetX()-0.5;
  this->_y2 = this->_input->GetY()-0.5;
  this->_z2 = this->_input->GetZ()-0.5;
}

double irtkShapeBasedInterpolateImageFunction::EvaluateInside(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  return this->_input->GetAsDouble(i, j, k, l);
}

double irtkShapeBasedInterpolateImageFunction::Evaluate(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  if ((i < 0) || (i >= this->_x) || (j < 0) || (j >= this->_y) || (k < 0) || (k >= this->_z)) {
    return this->_DefaultValue;
  } else {
    return this->_input->GetAsDouble(i, j, k, l);
  }
}

#endif

