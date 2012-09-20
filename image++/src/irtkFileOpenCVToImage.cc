/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkFileToImage.h>

#ifdef HAS_OPENCV

#include <irtkImageToOpenCv.h>

const char *irtkFileOpenCVToImage::NameOfClass()
{
  return "irtkFileOpenCVToImage";
}

int irtkFileOpenCVToImage::CheckHeader(const char *filename)
{
    IplImage *pimage = NULL;
    pimage= cvLoadImage(filename,0);
    if(pimage != NULL){
        cvReleaseImage(&pimage);
        return 1;
    }else{
        return 0;
    }
}

void irtkFileOpenCVToImage::SetInput(const char *filename){
    _pimage = NULL;
    _pimage= cvLoadImage(filename,0);

    // Read header
    this->ReadHeader();
}

irtkImage * irtkFileOpenCVToImage::GetOutput(){
    irtkGreyImage *output;
    irtkImageToOpenCv<irtkGreyPixel> itocg;
    itocg.SetOutput(_pimage);
    itocg.Invert();
    output = itocg.GetInput();
    return (irtkImage*)output;
}

void irtkFileOpenCVToImage::ReadHeader()
{
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
}

#endif

