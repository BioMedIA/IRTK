/*=========================================================================
 
 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id$
 Copyright : Imperial College, Department of Computing
 Visual Information Processing (VIP), 2008 onwards
 Date      : $Date$
 Version   : $Revision$
 Changes   : $Author$
 
 =========================================================================*/

#ifdef HAS_NIFTI

#include <irtkImage.h>

#include <irtkImageToFile.h>

#include <irtkNIFTI.h>

irtkImageToFileNIFTI::irtkImageToFileNIFTI() : irtkImageToFile()
{
	this->_headername = NULL;
}

irtkImageToFileNIFTI::~irtkImageToFileNIFTI()
{
	if (this->_headername != NULL) free(this->_headername);
}

void irtkImageToFileNIFTI::SetOutput(const char *name)
{
	int length;
	char *imagename;

	if (this->_headername != NULL) free(this->_headername);
	if (strstr(name, ".gz") == NULL) {
		this->_headername = strdup(name);
		imagename = strdup(name);
		length = strlen(name);
		imagename[length-1] = 'i';
		imagename[length-2] = 'i';
		imagename[length-3] = 'n';
	} else {
		this->_headername = strdup(name);
		imagename = strdup(name);
		length = strlen(name);
		imagename[length-4] = 'i';
		imagename[length-5] = 'i';
		imagename[length-6] = 'n';
	}
	this->irtkImageToFile::SetOutput(imagename);
}

const char *irtkImageToFileNIFTI::NameOfClass()
{
	return "irtkImageToFileNIFTI";
}

void irtkImageToFileNIFTI::Finalize()
{
	// Empty method
	return;
}

void irtkImageToFileNIFTI::Initialize()
{
  int x, y, z, t;
  double xsize, ysize, zsize, tsize;
  irtkMatrix i2w;
  double torigin;

  // Get image info
  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  t = this->_input->GetT();

  this->_input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);

  i2w = this->_input->GetImageToWorldMatrix();
  torigin = this->_input->ImageToTime(0);

  // Init header
  switch (this->_input->GetScalarType()) {
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_UINT8, i2w, torigin);
        break;
      }
    case IRTK_VOXEL_SHORT: {
        _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_INT16, i2w, torigin);
        break;
      }
    case IRTK_VOXEL_UNSIGNED_SHORT: {
        _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_UINT16, i2w, torigin);
        break;
      }
    case IRTK_VOXEL_INT: {
        _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_INT32, i2w, torigin);
        break;
      }
    case IRTK_VOXEL_UNSIGNED_INT: {
        _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_UINT32, i2w, torigin);
        break;
      }
    case IRTK_VOXEL_FLOAT: {
        _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_FLOAT32, i2w, torigin);
        break;
      }
    case IRTK_VOXEL_DOUBLE: {
        _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_FLOAT64, i2w, torigin);
        break;
      }
    default:
      cerr << "irtkImageToFileNIFTI::Run(): Unknown voxel type" << endl;
      exit(1);
  }
}

void irtkImageToFileNIFTI::Run()
{
	// Initialize filter
	this->Initialize();

	// Set filename in _hdr
	struct nifti_1_header nhdr = nifti_convert_nim2nhdr(_hdr.nim);
	_hdr.nim = nifti_convert_nhdr2nim(nhdr, this->_output);// This sets fname and iname
	_hdr.nim->iname_offset = 352;// Some nifti versions lose this on the way!

	/// Set data pointer in nifti image struct
	_hdr.nim->data = this->_input->GetScalarPointer();

	// Write hdr and data
	nifti_image_write(_hdr.nim);

	// Finalize filter
	this->Finalize();
}

#endif
