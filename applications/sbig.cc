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
#include <irtkResampling.h>
#include <irtkEuclideanDistanceTransform.h>
#include <irtkImageFunction.h>

void distancemap(irtkGreyImage &input, irtkRealImage &outputA)
{

  int x, y, z;
  irtkRealImage inputA, inputB, outputB;

  irtkEuclideanDistanceTransform<irtkRealPixel> *edt = NULL;
  edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
        (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform2D);

  inputA = input;
  inputB = input;

  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
        if (input(x, y, z) > 0.5) {
          inputA(x, y, z) = 1;
          inputB(x, y, z) = 0;
        } else {
          inputA(x, y, z) = 0;
          inputB(x, y, z) = 1;
        }
      }
    }
  }

  // Calculate EDT
  // cout << "Doing outside DT" << endl;
  edt->SetInput (& inputA);
  edt->SetOutput(&outputA);
  edt->Run();
  // cout << "Doing inside DT" << endl;
  edt->SetInput (& inputB);
  edt->SetOutput(&outputB);
  edt->Run();

  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
        outputA(x, y, z)  = sqrt(outputA(x, y, z)) - sqrt(outputB(x, y, z));
      }
    }
  }

  // Write image
  //  outputA.Write(output_name);

  // return 0;
}


irtkGreyImage threshold(irtkRealImage image, float minthres, float maxthres, short grey)
{

  int x, y, z;
  irtkGreyImage out;
  out = image;

  for (z = 0; z < image.GetZ(); z++) {
    for (y = 0; y < image.GetY(); y++) {
      for (x = 0; x < image.GetX(); x++) {
        if (image(x, y, z) <= maxthres&&image(x, y, z)>=minthres&&image(x, y, z)!=0)
          out(x, y, z) = grey;
        else out(x, y, z) = 0;
      }
    }
  }

  return(out);
}

irtkGreyImage threshold(irtkGreyImage image, short minthres, short maxthres, short grey)
{

  int x, y, z;

  for (z = 0; z < image.GetZ(); z++) {
    for (y = 0; y < image.GetY(); y++) {
      for (x = 0; x < image.GetX(); x++) {
        if (image(x, y, z) <= maxthres&&image(x, y, z)>=minthres)
          image(x, y, z) = grey;
        else
          image(x, y, z) = 0;
      }
    }
  }

  return(image);
}

irtkGreyImage AddImage(irtkGreyImage image1, irtkGreyImage image2, irtkGreyImage original_image)
{
  int x, y, z;
  short mingrey, maxgrey;
  irtkGreyImage image, tmpimage;
  image = image1;

  for (z = 0; z < image1.GetZ(); z++) {
    tmpimage = image1.GetRegion(z, 0);
    tmpimage.GetMinMax(&mingrey, &maxgrey);
    for (y = 0; y < image1.GetY(); y++) {
      for (x = 0; x < image1.GetX(); x++) {
        if (image1(x, y, z)> 0 && image2(x, y, z)==0 && mingrey!=maxgrey)
          image(x, y, z)=image1(x, y, z);
        else if (image1(x,y,z)==0 && image2(x,y,z)>0 && mingrey!=maxgrey)
          image(x, y, z) = image2(x, y, z);
        else if (image1(x,y,z)>0 && image2(x,y,z)>0 && mingrey!=maxgrey) {
          if (image1(x,y,z)>image2(x,y,z))
            image(x, y, z)= image1(x,y,z);
          else
            image(x, y, z)= image2(x,y,z);
        } else
          image(x, y, z) = 0;
      }
    }
  }
  return(image);
}

irtkRealImage sim3Ddmap(irtkRealImage image, float slices)
{
  int x, y, z;
  irtkRealPixel m;
  irtkRealImage output, output2;
  output  = image;
  output2 = image;

  for (z = 0; z < image.GetZ()-1; z++) {
    m=-1000;

    {
      for (y = 0; y < image.GetY(); y++) {
        for (x = 0; x < image.GetX(); x++) {
          output(x, y, z+1) = ((-image(x, y, z)))-slices;
          if ((-image(x, y, z+1)) < 0)
            m=max(m, output(x, y, z+1));
        }
      }
      for (y = 0; y < output.GetY(); y++) {
        for (x = 0; x < output.GetX(); x++) {
          output(x, y, z+1)=(output(x, y, z+1)-m)*slices/abs(slices+m);
          if (((-image(x, y, z+1)))<=0)
            output(x, y, z+1)=max((-image(x, y, z+1)), output(x, y, z+1));
          else
            output(x, y, z+1)=min((-image(x, y, z+1)), output(x, y, z+1));
        }
      }

    }

    {
      m=-1000;
      for (y = 0; y < image.GetY(); y++) {
        for (x = 0; x < image.GetX(); x++) {
          output(x, y, z) = (-image(x, y, z+1))-slices;
          if ((-image(x, y, z)) < 0)
            m=max(m, output(x, y, z));
        }
      }
      for (y = 0; y < output.GetY(); y++) {
        for (x = 0; x < output.GetX(); x++) {
          output(x, y, z)=(output(x, y, z)-m)*slices/abs(slices+m);
          if ((-image(x, y, z))<=0)
            output(x, y, z)=max((-image(x, y, z)), output(x, y, z));
          else
            output(x, y, z+1)=min((-image(x, y, z)), output(x, y, z));
        }
      }

    }
  }

  for (z = 0; z < output.GetZ(); z++) {
    for (y = 0; y < output.GetY(); y++) {
      for (x = 0; x < output.GetX(); x++) {
        output2(x,y,z) = -output(x,y,z);
      }
    }
  }

  return(output2);

}
/******************************************************************/


int main(int argc, char *argv[])
{
  int x, y, z;
  irtkGreyPixel min, max;
  float slices;
  short mingrey, maxgrey, grey; //min, max;
  irtkImageAttributes attr;
  irtkGreyImage imageA, tmpimageA, tmpimageB, tmpimage, totimage;
  irtkRealImage imageB, imageC, dmapimageB, imageoutputB, imageoutputC;

  if (argc != 4) {
    cerr << "Usage: sbig [input] [output] [slices]" << endl;
    exit(1);
  }

  imageA.Read(argv[1]);
  slices = atof(argv[3]);
  attr = imageA.GetImageAttributes();
  imageB.Initialize(attr);
  imageC.Initialize(attr);

  irtkImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
  irtkResampling<irtkRealPixel> resampling(attr._dx, attr._dy, (attr._dz*attr._z)/slices);

  imageA.GetMinMax(&mingrey, &maxgrey);
  cout << mingrey << " " << maxgrey << endl;

  for (grey=maxgrey; grey>0; grey--) {
    cout << grey << endl;
    tmpimageA = threshold(imageA, grey, grey, 1);
    tmpimageA.GetMinMax(&min, &max);
    if (min!=max)
      distancemap(tmpimageA, imageB);
    imageC=sim3Ddmap(imageB, slices/attr._z); //does not work yet

    for (z = 0; z < imageB.GetZ(); z++) {
      for (y = 0; y < imageB.GetY(); y++) {
        for (x = 0; x < imageB.GetX(); x++) {
          if (imageC(x,y,z) < imageB(x,y,z))
            imageB(x,y,z) = imageC(x,y,z);
        }
      }
    }

    resampling.SetInput(&imageB);
    resampling.SetOutput(&imageoutputB);
    resampling.SetInterpolator(interpolator);
    resampling.Run();

    tmpimage = threshold(imageoutputB, -500, 01, grey);
    if (totimage.GetX()==1) {
      totimage = tmpimage;
    } else {
      totimage = AddImage(tmpimage, totimage, imageA);
    }
    totimage.Write(argv[2]);
  }

  return 0;
}




