/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageToOpenCv.cc 101 2009-11-04 16:27:54Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-11-04 16:27:54 +0000 (‰∏? 04 ÂçÅ‰∏ÄÊú?2009) $
  Version   : $Revision: 101 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>

#ifdef HAS_OPENCV

template <class VoxelType> irtkImageToOpenCv<VoxelType>::irtkImageToOpenCv()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;
  min = 0;
  max = 0;
}

template <class VoxelType> irtkImageToOpenCv<VoxelType>::~irtkImageToOpenCv()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;
  min = 0;
  max = 0;
}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageToOpenCv::SetInput: Input is not an image\n";
    exit(1);
  }
  _input->GetMinMax(&min,&max);
}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::SetOutput(IplImage *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkImageToOpenCv::SetOutput: Output is not an opencv image\n";
    exit(1);
  }
}

template <class VoxelType> irtkGenericImage<VoxelType>* irtkImageToOpenCv<VoxelType>::GetInput(){
	return _input;
}

template <class VoxelType> IplImage* irtkImageToOpenCv<VoxelType>::GetOutput(){
	return _output;
}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::Initialize()
{
  // Check inputs and outputs
  if (_input == NULL && _output == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input nor output" << endl;
    exit(1);
  }

   if (_input == NULL) {
   GenInput();
  }

  if (_output == NULL) {
   GenOutput();
  }

}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::Finalize()
{
  
}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::GenOutput()
{
	if(!_output)
		_output= cvCreateImage(cvSize(_input->GetX(), _input->GetY()), IPL_DEPTH_8U, 1);
	else
		cerr<<this->NameOfClass() << "::GenOutput: Output is not empty" << endl;
}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::GenInput()
{
	if(!_input){
		irtkImageAttributes attr;
		attr._yaxis[0] = 0; attr._yaxis[1] = -1; attr._yaxis[2] = 0;
		attr._xaxis[0] = -1; attr._xaxis[1] = 0; attr._xaxis[2] = 0;
		attr._x = _output->width;
		attr._y = _output->height;
		_input= new irtkGenericImage<VoxelType>(attr);
	}else
		cerr<<this->NameOfClass() << "::GenOutput: Input is not empty" << endl;
}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::Run(int k)
{
  int i,j,tmp;

  // Do the initial set up
  this->Initialize();

  for (j = 0; j < _input->GetY(); j++) {
	  for (i = 0; i < _input->GetX(); i++) {
		  if(max > 255){
		  tmp = (_input->GetAsDouble(i,j,k,0) * 255 /max);
		  }else{
		  tmp = (_input->GetAsDouble(i,j,k,0));
		  }
		  //*cvPtr2D(_output,j,i) = tmp;
		  _output->imageData[j*_output->widthStep + i] = tmp;
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

template <class VoxelType> void irtkImageToOpenCv<VoxelType>::Invert(int k)
{
  int i,j;
  int tmp;

  // Do the initial set up
  this->Initialize();

  for (j = 0; j < _input->GetY(); j++) {
	  for (i = 0; i < _input->GetX(); i++) {
		  tmp = *cvPtr2D(_output,j,i);
		  _input->PutAsDouble(i,j,k,tmp);
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkImageToOpenCv<unsigned char>;
template class irtkImageToOpenCv<short>;
template class irtkImageToOpenCv<unsigned short>;
template class irtkImageToOpenCv<float>;
template class irtkImageToOpenCv<double>;

#endif