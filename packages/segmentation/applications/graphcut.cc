/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkSegmentationFunction.h>

int Zcheck = 0;
char *outputname = NULL;

void usage()
{
  cerr << "Usage: graphcut [NumberOfImages] [Input_1...Input_n]" << endl;
  cerr << "-GMM [value]                       GMM number for EM"<<endl;
  cerr << "-atlas [n atlas1...atlasn]         atlas information for segmentation, no need to provide background" << endl;
  cerr << "-numberofcomponents [n_1...n_n+1]  number of components per atlas + number of components of background must be used with atlas."<<endl;
  cerr << "-graphcutmode [value]              Graphcut mode 1D 2D 3D 4D (1,2,3,4 with 4 connective neighbors)"<<endl;
  cerr << "-connect 0/1	                      0 neighbor connection 1 cubic connection"<<endl;
  cerr << "-dataweight [value]                Data term weight for graphcut"<<endl;
  cerr << "-timeweight [value]                Time weight for boundary term in 4D graphcut case" <<endl;
  cerr << "-outputname [name]                 Segmentation output name" << endl;
  exit(1);
}

int main( int argc, char** argv )
{
  irtkRealImage late;
  char buffer[255];
  irtkRealImage **atlas = NULL;
  irtkGreyImage sub,output;
  irtkImageAttributes atr;
  irtkSegmentationFunction cf;
  int i,j,k,l,ok,*c = NULL;
  int numberOfImages = 0, sublate = 0, numberofatlas = 0, cutmode = 3, seperate = 0, connect = 0;
  double dataweight = 0.5,timeweight = 1.0;
  int gmm = 12;

  if( argc < 3 ) usage();

  // Number Of Images
  numberOfImages = atoi(argv[1]);
  argc--;
  argv++;
  irtkGreyImage **input = new irtkGreyImage *[numberOfImages];

  // Read image
  for (i = 0; i < numberOfImages; i++) {
    cout << argv[1] << endl;
    input[i] = new irtkGreyImage(argv[1]);
    argv++;
    argc--;
  }

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-sub") == 0)) {
      argc--;
      argv++;
      late.Read(argv[1]);
      sublate = 1;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-atlas") == 0)) {
      argc--;
      argv++;
      numberofatlas = atoi(argv[1]);
      atlas = new irtkRealImage*[numberofatlas];
      argc--;
      argv++;
      // Read atlas for each tissue
      for (i = 0; i < numberofatlas; i++) {
        atlas[i] = new irtkRealImage;
        atlas[i]->Read(argv[1]);
        cerr << "Image " << i <<" = " << argv[1] <<endl;
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-numberofcomponents") == 0)) {
      argc--;
      argv++;
      c = new int[numberofatlas+1];
      // Read atlas for each tissue
      for (i = 0; i < numberofatlas+1; i++) {
        c[i] = atoi(argv[1]);
        cerr << "Component " << i <<" has " << argv[1] <<" GMM models" <<endl;
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-outputname") == 0)) {
      argc--;
      argv++;
      outputname = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-connect") == 0)) {
      argc--;
      argv++;
      connect = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-timeweight") == 0)) {
      argc--;
      argv++;
      timeweight = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-graphcutmode") == 0)) {
      argc--;
      argv++;
      cutmode = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-seperate") == 0)) {
      argc--;
      argv++;
      seperate = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dataweight") == 0)) {
      argc--;
      argv++;
      dataweight = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-GMM") == 0)) {
      argc--;
      argv++;
      gmm = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  //find out number of time frame
  atr = input[0]->GetImageAttributes();

  //generate integral image
  sub.Initialize(atr);

  // generate output file name
  if(!outputname)
    sprintf(buffer,"test.nii");
  else
    sprintf(buffer,"%s",outputname);

  output.Initialize(atr);

  //get threshold
  if(sublate) {
    for(l = 0;l<atr._t;l++) {
      for(k = 0; k < atr._z; k++) {
        for(j = 0; j < atr._y; j++) {
          for(i = 0; i < atr._x; i++) {
            sub.PutAsDouble(i,j,k,l,input[0]->GetAsDouble(i,j,k,l));
          }
        }
      }
    }
    for(l = 0;l<atr._t;l++) {
      for(k = 0; k < atr._z; k++) {
        for(j = 0; j < atr._y; j++) {
          for(i = 0; i < atr._x; i++) {
            int temp;
            temp = sub.GetAsDouble(i,j,k,l) - late.GetAsDouble(i,j,k,l);
            if (temp < sub.GetAsDouble(i,j,k,l)*2.0/3.0) temp = sub.GetAsDouble(i,j,k,l)*2.0/3.0;
            input[0]->PutAsDouble(i,j,k,l,temp);
          }
        }
      }
    }
  }
  //input[0]->Write(buffer);
  //blue before segmentation
  cout << "Blurring for segmentation ... "; cout.flush();
  irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(1, 0);
  blurring.SetInput (input[0]);
  blurring.SetOutput(input[0]);
  blurring.Run();
  cout << "done" << endl;
  //evaluate threshold
  if(numberofatlas) {
    cf.EvaluateGraphCut(&output, input,atlas,numberOfImages,numberofatlas, timeweight,cutmode,connect,dataweight,0,c);
  } else {
    cf.EvaluateGraphCut(&output, input,numberOfImages,timeweight,cutmode,connect,dataweight,0,gmm);
  }

  output.Write(buffer);

  for (i = 0; i < numberOfImages; i++) {
    delete input[i];
  }
  delete input;
  if(numberofatlas > 0) {
    delete []atlas;
  }
  if(c) {
    delete []c;
  }
}
