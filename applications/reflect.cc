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

char *input_name = NULL, *output_name = NULL;

char **reflect_list = NULL;

void usage()
{
  cerr << "Usage: reflect [in] [out] [reflection_1] <reflection_2> .. <reflection_n>" << endl;
  cerr << "Where the reflections are chosen from:" << endl;
  cerr << " " << endl;
  cerr << "        -x -y -z -xy -xz -yz" << endl;
  cerr << " " << endl;
  cerr << "The reflections are processed in the order given." << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, count;
  irtkBaseImage *image;

  if (argc < 4){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  count = argc - 1;
  reflect_list = new char *[count];
  for (i = 0; i < count; ++i){
    reflect_list[i] = argv[1];
    argc--;
    argv++;
  }

  // Read input
  irtkFileToImage *reader = irtkFileToImage::New(input_name);
  image = reader->GetOutput();

  for (i = 0; i < count; ++i){
    if (strcmp("-x", reflect_list[i]) == 0){
      image->ReflectX();
    }

    if (strcmp("-y", reflect_list[i]) == 0){
     image->ReflectY();
    }

    if (strcmp("-z", reflect_list[i]) == 0){
        image->ReflectZ();
    }

    if (strcmp("-xy", reflect_list[i]) == 0 ||
        strcmp("-yx", reflect_list[i]) == 0){
      image->FlipXY();
    }

    if (strcmp("-xz", reflect_list[i]) == 0 ||
        strcmp("-zx", reflect_list[i]) == 0){
      image->FlipXZ();
    }

    if (strcmp("-yz", reflect_list[i]) == 0 ||
        strcmp("-zy", reflect_list[i]) == 0){
      image->FlipYZ();
    }
  }

  // Write region
  image->Write(output_name);

  delete [] reflect_list;

  return 0;
}
