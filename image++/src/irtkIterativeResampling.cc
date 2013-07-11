/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkIterativeResampling.h>

template <class VoxelType> irtkIterativeResampling<VoxelType>::irtkIterativeResampling(double new_xsize, double new_ysize, double new_zsize) : irtkResampling<VoxelType>(new_xsize, new_ysize, new_zsize)
{}

template <class VoxelType> const char *irtkIterativeResampling<VoxelType>::NameOfClass()
{
  return "irtkIterativeResampling";
}

template <class VoxelType> void irtkIterativeResampling<VoxelType>::Run()
{
  this->irtkResampling<VoxelType>::Initialize();

  if (this->_input->GetZ() == 1) {
	this->Run2D();
  } else {
	this->Run3D();
  }

  // Do the final cleaning up
  this->Finalize();
}

template <class VoxelType> void irtkIterativeResampling<VoxelType>::Run2D()
{
  double xsize, ysize, zsize;
  irtkGenericImage<VoxelType> tmp_in, tmp_out;

  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  tmp_in = *this->_input;
  tmp_out = *this->_input;

  while ((xsize < this->_XSize) && (ysize < this->_YSize)) {
	xsize *= 2.0;
	ysize *= 2.0;

    _innerResample = new irtkResampling<VoxelType>(xsize, ysize, this->_ZSize);
	_innerResample->SetInput(&tmp_in);
	_innerResample->SetOutput(&tmp_out);
	_innerResample->SetInterpolator(this->_Interpolator);
	_innerResample->Run();

	//Copy output as new input for next iteration
	tmp_in = tmp_out;

	delete _innerResample;
  }

  //Do final resampling if necessary to have everything OK
  if ((xsize != this->_XSize) || (ysize != this->_YSize)) {
    _innerResample = new irtkResampling<VoxelType>(this->_XSize, this->_YSize, this->_ZSize);
    _innerResample->SetInput(&tmp_in);
    _innerResample->SetOutput(&tmp_out);
    _innerResample->SetInterpolator(this->_Interpolator);
    _innerResample->Run();

    delete _innerResample;
  }

  //Copy into output
  *this->_output = tmp_out;
}

template <class VoxelType> void irtkIterativeResampling<VoxelType>::Run3D()
{
  double xsize, ysize, zsize;
  irtkGenericImage<VoxelType> tmp_in, tmp_out;

  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  tmp_in = *this->_input;
  tmp_out = *this->_input;

  while ((xsize < this->_XSize) && (ysize < this->_YSize) && (zsize < this->_ZSize)) {
	xsize *= 2.0;
	ysize *= 2.0;
	zsize *= 2.0;

    _innerResample = new irtkResampling<VoxelType>(xsize, ysize, zsize);
	_innerResample->SetInput(&tmp_in);
	_innerResample->SetOutput(&tmp_out);
	_innerResample->SetInterpolator(this->_Interpolator);
	_innerResample->Run();

	//Copy output as new input for next iteration
	tmp_in = tmp_out;

	delete _innerResample;
  }

  //Do final resampling if necessary to have everything OK
  if ((xsize != this->_XSize) || (ysize != this->_YSize) || (zsize != this->_ZSize)) {
    _innerResample = new irtkResampling<VoxelType>(this->_XSize, this->_YSize, this->_ZSize);
    _innerResample->SetInput(&tmp_in);
    _innerResample->SetOutput(&tmp_out);
    _innerResample->SetInterpolator(this->_Interpolator);
    _innerResample->Run();

    delete _innerResample;
  }

  //Copy into output
  *this->_output = tmp_out;
}

template class irtkIterativeResampling<char>;
template class irtkIterativeResampling<unsigned char>;
template class irtkIterativeResampling<short>;
template class irtkIterativeResampling<unsigned short>;
template class irtkIterativeResampling<int>;
template class irtkIterativeResampling<float>;
template class irtkIterativeResampling<double>;
