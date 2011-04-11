#include <irtkSegmentationFunction.h>

int Zcheck = 0;
char *outputname = NULL;

void usage()
{
  cerr << "Usage: graphcut [NumberOfImages] [Input_1...Input_n]" << endl;
  cerr << "-GMM [value,value]                 GMM num for EM,num of GMM Region in Background of EM"<<endl;
  cerr << "-atlas [n atlas1...atlasn]         atlas information for segmentation, the first one is the foreground atlas, rest is back ground atlas" << endl;
  cerr << "-numberofcomponents [n1...nn]      number of components per atlas can't be used when do not have atlas."<<endl;
  cerr << "-graphcutmode [value]              Graphcut mode 1D 2D 3D 4D (1,2,3,4 with 4 connective neighbors)"<<endl;
  cerr << "-connect 0/1	                      0 neighbor connection 1 cubic connection"<<endl;
  cerr << "-regionweight [value]              Region term weight for graphcut"<<endl;
  cerr << "-timeweight [value]                Time weight for boundary term in 4D graphcut case" <<endl;
  cerr << "-seperate                          Seperate atlas and run seperate ems segmentation on all time frames"<<endl;
  cerr << "-outputname [name]                 Segmentation output name" << endl;
  cerr << "-outputlabel [value]               Output Label in the output segmentation default 1" << endl;
  cerr << "-load                              Load segmentation result" << endl;
  exit(1);
}

int main( int argc, char** argv )
{
  irtkRealImage late;
  char buffer[255];
  irtkRealImage **atlas = NULL;
  irtkRealImage threshold,sub,output;
  irtkImageAttributes atr;
  irtkSegmentationFunction cf;
  int i,j,k,l,ok,*c = NULL;
  int numberOfImages = 0, sublate = 0, numberofatlas = 0, cutmode = 3, seperate = 0, connect = 0, outputvalue = 1,loadsegmentation = 0;
  double regionweight = 0.5,timeweight = 1.0;
  int gmm = 12, bgmm = 4;

  if( argc < 3 ) usage();

  // Number Of Images
  numberOfImages = atoi(argv[1]);
  argc--;
  argv++;
  irtkRealImage **input = new irtkRealImage *[numberOfImages];

  // Read image
  for (i = 0; i < numberOfImages; i++) {
    cout << argv[1] << endl;
    input[i] = new irtkRealImage(argv[1]);
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
    if ((ok == false) && (strcmp(argv[1], "-load") == 0)) {
      argc--;
      argv++;
      loadsegmentation = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-regionweight") == 0)) {
      argc--;
      argv++;
      regionweight = atof(argv[1]);
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
      bgmm = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-outputlabel") == 0)) {
      argc--;
      argv++;
      outputvalue = atoi(argv[1]);
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

  if(loadsegmentation)
    output.Read(buffer);
  else
    output.Initialize(atr);

  threshold.Initialize(atr);

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
  irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(1, 0);
  blurring.SetInput (input[0]);
  blurring.SetOutput(input[0]);
  blurring.Run();
  cout << "done" << endl;
  //evaluate threshold
  if(numberofatlas) {
    cf.EvaluateGraphCut(&threshold, input,atlas,numberOfImages,numberofatlas, timeweight,cutmode,connect,regionweight,seperate,0,c);
  } else {
    cf.EvaluateGraphCut(&threshold, input,numberOfImages,timeweight,cutmode,connect,regionweight,0,gmm,bgmm);
  }

  for(int l = 0; l < threshold.GetT(); l++) {
    for(int k = 0; k < threshold.GetZ(); k++) {
      for(int j = 0; j < threshold.GetY(); j++) {
        for(int i = 0; i < threshold.GetX(); i++) {
          if(threshold.GetAsDouble(i,j,k,l) > 0) {
            output.PutAsDouble(i,j,k,l,outputvalue);
          }
        }
      }
    }
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
