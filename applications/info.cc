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

char *input_name = NULL;

void usage()
{
  cerr << "Usage: info [image] [-real] [-grey] [-byte]\n";
  exit(1);
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    usage();
  }

  if (argc == 3) {
    if (strcmp(argv[2], "-byte") == 0) {
      irtkBytePixel min, max;
      irtkByteImage image;
      irtkFileToImage<irtkBytePixel> *reader =
        irtkFileToImage<irtkBytePixel>::New(argv[1]);
      cout << "Information from irtkFileToImage" << endl;
      reader->Print();
      cout << "Information from irtkGenericImage" << endl;
      image.Read(argv[1]);
      image.Print();
      image.GetMinMax(&min, &max);
      cout << min << " " << max << endl;
      cout << "Image to world matrix" << endl;
      image.GetImageToWorldMatrix().Print();
      cout << "World to image matrix" << endl;
      image.GetWorldToImageMatrix().Print();
    } else {
      if (strcmp(argv[2], "-grey") == 0) {
        irtkGreyPixel min, max;
        irtkGreyImage image;
        irtkFileToImage<irtkGreyPixel> *reader =
          irtkFileToImage<irtkGreyPixel>::New(argv[1]);
        cout << "Information from irtkFileToImage" << endl;
        reader->Print();
        cout << "Information from irtkGenericImage" << endl;
        image.Read(argv[1]);
        image.Print();
        image.GetMinMax(&min, &max);
        cout << min << " " << max << endl;
        cout << "Image to world matrix" << endl;
        image.GetImageToWorldMatrix().Print();
        cout << "World to image matrix" << endl;
        image.GetWorldToImageMatrix().Print();
      } else {
        if (strcmp(argv[2], "-real") == 0) {
          irtkRealPixel min, max;
          irtkRealImage image;
          irtkFileToImage<irtkRealPixel> *reader =
            irtkFileToImage<irtkRealPixel>::New(argv[1]);
          cout << "Information from irtkFileToImage" << endl;
          reader->Print();
          cout << "Information from irtkGenericImage" << endl;
          image.Read(argv[1]);
          image.Print();
          image.GetMinMax(&min, &max);
          cout << min << " " << max << endl;
          cout << "Image to world matrix" << endl;
          image.GetImageToWorldMatrix().Print();
          cout << "World to image matrix" << endl;
          image.GetWorldToImageMatrix().Print();
        } else {
          usage();
        }
      }
    }
  } else {
    irtkGreyPixel min, max;
    irtkGreyImage image;
    irtkFileToImage<irtkGreyPixel> *reader =
      irtkFileToImage<irtkGreyPixel>::New(argv[1]);
    cout << "Information from irtkFileToImage" << endl;
    reader->Print();
    cout << "Information from irtkGenericImage" << endl;
    image.Read(argv[1]);
    image.Print();
    image.GetMinMax(&min, &max);
    cout << min << " " << max << endl;
    cout << "Image to world matrix" << endl;
    image.GetImageToWorldMatrix().Print();
    cout << "World to image matrix" << endl;
    image.GetWorldToImageMatrix().Print();
  }

  return 0;
}
