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

#include <irtkImageToFile.h>

#ifdef HAS_PNG

#include <png.h>

template <class VoxelType> const char *irtkImageToFilePNG<VoxelType>::NameOfClass()
{
  return "irtkImageToFilePNG";
}

template <class VoxelType> void irtkImageToFilePNG<VoxelType>::Initialize()
{
  // Initialize base class
  this->irtkImageToFile<VoxelType>::Initialize();

  if (this->_input->GetZ() != 1) {
    cerr << this->NameOfClass() << " supports only 2D images" << endl;
    exit(1);
  }
}

template <> void irtkImageToFilePNG<irtkRGBPixel>::Run()
{
  int x, y;

  // Initialize filter
  this->Initialize();

  png_structp png_ptr = png_create_write_struct
                        (PNG_LIBPNG_VER_STRING, (png_voidp)NULL, NULL, NULL);

  if (!png_ptr) {
    cerr << "irtkImageToFilePNG<irtkRGBPixel>::Run: Unable to write PNG file!" << endl;
    exit(1);
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr,
                             (png_infopp)NULL);
    cerr << "irtkImageToFilePNG<irtkRGBPixel>::Run: Unable to write PNG file!" << endl;;
    exit(1);
  }

  if (setjmp(png_ptr->jmpbuf)) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    exit(1);
  }

  // Open file
  FILE *fp = fopen(_output, "wb");

  // Initialize PNG I/O
  png_init_io(png_ptr, fp);

  // Initialize header
  png_set_IHDR(png_ptr, info_ptr, _input->GetX(), _input->GetY(),
               8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);

  // Write header
  png_write_info(png_ptr, info_ptr);

  // Copy image
  png_byte *data = new png_byte  [3*_input->GetX()*_input->GetY()];
  png_byte **ptr = new png_byte *[_input->GetY()];

  png_byte *ptr2data = data;
  for (y = 0; y < _input->GetY(); y++) {
    for (x = 0; x < _input->GetX(); x++) {
      *ptr2data = _input->Get(x, y, 0, 0)._r;
      ptr2data++;
      *ptr2data = _input->Get(x, y, 0, 0)._g;
      ptr2data++;
      *ptr2data = _input->Get(x, y, 0, 0)._b;
      ptr2data++;
    }
    ptr[_input->GetY() - y - 1] = &(data[_input->GetX()*y*3]);
  }
  png_write_image(png_ptr, ptr);
  png_write_end(png_ptr, info_ptr);

  // Delete pointers
  delete [] data;
  delete [] ptr;

  // Destroy PNG data
  png_destroy_write_struct(&png_ptr, &info_ptr);

  // Close file
  fclose(fp);

  // Finalize filter
  this->Finalize();
}

template class irtkImageToFilePNG<irtkRGBPixel>;

#endif
