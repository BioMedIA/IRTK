/*=========================================================================

  Date      : $Date: 29.04.2010$
  Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$

=========================================================================*/

#include <irtkLargeDeformationSciCalcPack.h>





/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                            THE TOOLS
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ image transport +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void transportImage(char * ImageInit,char * ImageFinal,int TimeInit,int TimeFinal,char * VFX,char * VFY,char * VFZ){
  ScalarField ImInit;
  ScalarField ImFinal;
  VectorField VelocityField;
  VectorField Map;
  VectorField IdMap;
  int x,y,z;
  
  //init
  ImInit.Read(ImageInit);
  ImFinal.Read(ImageInit);  //to allocate it
  VelocityField.Read(VFX,VFY,VFZ);
  IdMap.CreateVoidField(VelocityField.NX,VelocityField.NY,VelocityField.NZ);
  
  for (z = 0; z < IdMap.NZ; z++)  for (y = 0; y < IdMap.NY; y++) for (x = 0; x < IdMap.NX; x++){
    IdMap.P(static_cast<float>(x),0,x,y,z);
    IdMap.P(static_cast<float>(y),1,x,y,z);
    IdMap.P(static_cast<float>(z),2,x,y,z);
  }

  
  //compute the mapping
  CptMappingFromVeloField(TimeInit,&IdMap,&VelocityField,&Map);
  
  //project the image
  Project3Dimage(&ImInit,&Map,&ImFinal,TimeFinal);
  
  //save the projected image
  ImFinal.Write(ImageFinal,ImageInit);

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ transport several images  +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//largedeformationsTools -MultiTransport [VFX][VFY][VFZ] [Src TimeSub] [Nb Ima] [Src Ima 1][Trg TimeSub 1][Trg Ima 1] ...
//-> Transport [Nb Ima] images [Src Ima i] from the time subdivision [Src TimeSub] to the time subdivision [Trg TimeSub i]
//through the velocity field [VFX][VFY][VFZ]. The transported images are saved in [Trg Ima i].
void MultiTransport(int argc, char **argv){
  char ImageInit[256];
  char ImageFinal[256];
  ScalarField ImInit;
  ScalarField ImFinal;
  VectorField VelocityField;
  VectorField Map;
  VectorField IdMap;
  int x,y,z;
  char VFX[256];
  char VFY[256];
  char VFZ[256];
  int NbImages;
  int TimeFinal;
  int TimeInit;
  int i;
  
  //read global parameters
  argc--; argv++;
  strcpy(VFX,argv[1]); //velocity field X
  argc--; argv++;
  strcpy(VFY,argv[1]); //velocity field Y
  argc--; argv++;
  strcpy(VFZ,argv[1]); //velocity field Z
  argc--; argv++;
  TimeInit = atoi(argv[1]);
  argc--;  argv++;
  NbImages = atoi(argv[1]);
  argc--;  argv++;
  
  //compute the mapping
  VelocityField.Read(VFX,VFY,VFZ);
  IdMap.CreateVoidField(VelocityField.NX,VelocityField.NY,VelocityField.NZ);
  
  for (z = 0; z < IdMap.NZ; z++)  for (y = 0; y < IdMap.NY; y++) for (x = 0; x < IdMap.NX; x++){
    IdMap.P(static_cast<float>(x),0,x,y,z);
    IdMap.P(static_cast<float>(y),1,x,y,z);
    IdMap.P(static_cast<float>(z),2,x,y,z);
  }

  CptMappingFromVeloField(TimeInit,&IdMap,&VelocityField,&Map);
  
  //transport all the images
  for (i=0;i<NbImages;i++){
    //read the parameters
    strcpy(ImageInit,argv[1]);
    argc--; argv++;
    TimeFinal = atoi(argv[1]);
    argc--;  argv++;
    strcpy(ImageFinal,argv[1]);
    argc--; argv++;
    
    //read and allocate the images
    ImInit.Read(ImageInit);
    ImFinal.Read(ImageInit);

    //project the image
    Project3Dimage(&ImInit,&Map,&ImFinal,TimeFinal);
    
    //save the projected image
    ImFinal.Write(ImageFinal,ImageInit);
  }
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ weighted sum +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int weightedSum(int argc, char **argv)
{
  ScalarField in;
  ScalarField out;
  int NbImages;
  double weight;
  char *input_name = NULL, *output_name = NULL;
  int x,y,z,t,i;
  
  
  argc--;
  argv++;
  NbImages = atoi(argv[1]);
  argc--;
  argv++;
  
  for (i=0;i<NbImages;i++){
    if (i==0){
      // first image ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);
      
      out.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      
      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.P(in.G(x,y,z,t)*weight,x,y,z,t);
    }
    else{
      // other images ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);
      
      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.Add(in.G(x,y,z,t)*weight,x,y,z,t);
    }
  }
  
  output_name = argv[1];
  argc--;
  argv++;
  out.Write(output_name,input_name);
  
  return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ weighted sum +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int weightedMean(int argc, char **argv)
{
  ScalarField in;
  ScalarField out;
  int NbImages;
  double weight;
  double SumWeight;
  char *input_name = NULL, *output_name = NULL;
  int x,y,z,t,i;
  
  
  argc--;
  argv++;
  NbImages = atoi(argv[1]);
  argc--;
  argv++;
  
  SumWeight=0.;
  
  for (i=0;i<NbImages;i++){
    if (i==0){
      // first image ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);
      
      out.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      
      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.P(in.G(x,y,z,t)*weight,x,y,z,t);
      
      SumWeight+=weight;
    }
    else{
      // other images ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);
      
      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.Add(in.G(x,y,z,t)*weight,x,y,z,t);
      
      SumWeight+=weight;
    }
  }
  
  //compute the mean
  for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
    out.P(out.G(x,y,z,t)/SumWeight,x,y,z,t);
  
  
  output_name = argv[1];
  argc--;
  argv++;
  out.Write(output_name,input_name);
  
  return 0;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++ weighted standard deviation ++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int weightedStd(int argc, char **argv)
{
  ScalarField in;
  ScalarField std;
  ScalarField mean;
  int NbImages;
  double weight;
  double SumWeight;
  char *input_name = NULL, *output_name = NULL;
  int x,y,z,t,i;
  
  argc--;
  argv++;
  NbImages = atoi(argv[1]);
  argc--;
  argv++;
  
  
  //1) compute the mean
  SumWeight=0.;
  
  for (i=0;i<NbImages;i++){
    //... read parameters
    weight = atof(argv[2*i+1]);
    input_name  = argv[2*i+2];
    
    //... do the treatment
    in.Read(input_name);
    
    if (i==0){
      mean.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      for (t = 0; t < mean.NT; t++) for (z = 0; z < mean.NZ; z++)  for (y = 0; y < mean.NY; y++) for (x = 0; x < mean.NX; x++)
        mean.P(in.G(x,y,z,t)*weight,x,y,z,t);
    }
    else{
      for (t = 0; t < mean.NT; t++) for (z = 0; z < mean.NZ; z++)  for (y = 0; y < mean.NY; y++) for (x = 0; x < mean.NX; x++)
        mean.Add(in.G(x,y,z,t)*weight,x,y,z,t);
    }
    
    SumWeight+=weight;
  }
  
  //divide by the sum of the weights to obtain a mean
  for (t = 0; t < mean.NT; t++) for (z = 0; z < mean.NZ; z++)  for (y = 0; y < mean.NY; y++) for (x = 0; x < mean.NX; x++)
    mean.P(mean.G(x,y,z,t)/SumWeight,x,y,z,t);
  
  
  //2) compute the standard deviation
  
  SumWeight=0.;
  
  for (i=0;i<NbImages;i++){
    //... read parameters
    weight = atof(argv[2*i+1]);
    input_name  = argv[2*i+2];
    
    //... do the treatment
    in.Read(input_name);
    
    if (i==0){
      std.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      for (t = 0; t < std.NT; t++) for (z = 0; z < std.NZ; z++)  for (y = 0; y < std.NY; y++) for (x = 0; x < std.NX; x++)
        std.P((in.G(x,y,z,t)-mean.G(x,y,z,t))*(in.G(x,y,z,t)-mean.G(x,y,z,t))*weight,x,y,z,t);
    }
    else{
      for (t = 0; t < std.NT; t++) for (z = 0; z < std.NZ; z++)  for (y = 0; y < std.NY; y++) for (x = 0; x < std.NX; x++)
        std.Add((in.G(x,y,z,t)-mean.G(x,y,z,t))*(in.G(x,y,z,t)-mean.G(x,y,z,t))*weight,x,y,z,t);
    }
    
    SumWeight+=weight;
  }
  
  //divide by the sum of the weights and compute the square root to obtain a std
  for (t = 0; t < std.NT; t++) for (z = 0; z < std.NZ; z++)  for (y = 0; y < std.NY; y++) for (x = 0; x < std.NX; x++)
    std.P(sqrt(std.G(x,y,z,t)/SumWeight),x,y,z,t);
  
  
  output_name = argv[2*NbImages+1];
  argc--;
  argv++;
  std.Write(output_name,input_name);
  
  return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++ sum of squared differences between two images +++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double MeasureSSD(int argc, char **argv)
{
  ScalarField image1;
  ScalarField image2;
  double SSD;
  int x,y,z,t;
  int nbPts;
  char File1[256];
  char File2[256];
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  
  //intialisation
  image1.Read(File1);
  image2.Read(File2);
  
  SSD=0.;
  nbPts=0;
  
  //computation
  for (t = 0; t < image1.NT; t++) for (z = 0; z < image1.NZ; z++)  for (y = 0; y < image1.NY; y++) for (x = 0; x < image1.NX; x++){
    SSD+=(image1.G(x,y,z,t)-image2.G(x,y,z,t))*(image1.G(x,y,z,t)-image2.G(x,y,z,t));
    nbPts++;
  }
  
  SSD=sqrt(SSD);
  SSD/=(double)nbPts;
  
  cout << "SSD=" << SSD << "\n";
  return SSD;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++ surfacic average and std values ++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double Surf_Av_Std_values(int argc, char **argv)
{
  ScalarField RefIma;
  ScalarField SalarsIma;
  double Thresh;
  int x,y,z,t;
  int nbPts;
  char File1[256];
  char File2[256];
  double av,std;
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  Thresh = atof(argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  
  //intialisation
  RefIma.Read(File1);
  SalarsIma.Read(File2);
  
  av=0.;
  nbPts=0;
  
  //computation
  for (t = 0; t < RefIma.NT; t++) for (z = 0; z < RefIma.NZ; z++)  for (y = 0; y < RefIma.NY; y++) for (x = 0; x < RefIma.NX; x++)
    if (RefIma.G(x,y,z,t)>Thresh)
      if ((RefIma.G(x+1,y,z,t)<Thresh)||(RefIma.G(x-1,y,z,t)<Thresh)||(RefIma.G(x,y+1,z,t)<Thresh)||(RefIma.G(x,y-1,z,t)<Thresh)||(RefIma.G(x,y,z+1,t)<Thresh)||(RefIma.G(x,y,z-1,t)<Thresh)){
        av+=SalarsIma.G(x,y,z,t);
        nbPts++;
      }
  
  av/=(double)nbPts;
  
  std=0;
  
  for (t = 0; t < RefIma.NT; t++) for (z = 0; z < RefIma.NZ; z++)  for (y = 0; y < RefIma.NY; y++) for (x = 0; x < RefIma.NX; x++)
    if (RefIma.G(x,y,z,t)>Thresh)
      if ((RefIma.G(x+1,y,z,t)<Thresh)||(RefIma.G(x-1,y,z,t)<Thresh)||(RefIma.G(x,y+1,z,t)<Thresh)||(RefIma.G(x,y-1,z,t)<Thresh)||(RefIma.G(x,y,z+1,t)<Thresh)||(RefIma.G(x,y,z-1,t)<Thresh)){
        std+=(SalarsIma.G(x,y,z,t)-av)*(SalarsIma.G(x,y,z,t)-av);
      }

  std/=(double)nbPts;
  std=sqrt(std);

  cout << "av=" << av <<  "   std=" << std << "\n";
  return av;
}






/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                        INPUTS MANAGEMENT
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void usage(){
  cerr << "Usage: largedeformationsTools <options>\n";
  cerr << "\n";
  cerr << "Where <options> are the following:\n";
  cerr << "\n";
  cerr << "-transport [Time subdivision i][Image T_i] [Velo Field X][Velo Field Y][Velo Field Z] [Time subdivision j][Image T_j]\n";
  cerr << "    -> Transport the image [Image T_i] from the time subdivision i to j using the velocity field. The result is\n";
  cerr << "    saved in [Image T_j].\n";
  cerr << "\n";
  cerr << "-MultiTransport [VFX][VFY][VFZ] [Src TimeSub] [Nb Ima] [Src Ima 1][Trg TimeSub 1][Trg Ima 1] ...\n";
  cerr << "    -> Transport [Nb Ima] images [Src Ima i] from the time subdivision [Src TimeSub] to the time subdivision [Trg TimeSub i]\n";
  cerr << "    through the velocity field [VFX][VFY][VFZ]. The transported images are saved in [Trg Ima i].\n";
  cerr << "\n";
  cerr << "-weightedSum [Nb Images] [Weight1][Image1] [Weight2][Image2]... [Output Image]\n";
  cerr << "    -> Compute the weighted sum of one or several images\n";
  cerr << "\n";
  cerr << "-weightedMean [Nb Images] [Weight1][Image1] [Weight2][Image2]... [Output Image]\n";
  cerr << "    -> Compute the weighted mean of several images\n";
  cerr << "\n";
  cerr << "-weightedStd [Nb Images] [Weight1][Image1] [Weight2][Image2]... [Output Image]\n";
  cerr << "    -> Compute the weighted standard deviation between several images\n";
  cerr << "\n";
  cerr << "-SSD [Image1] [Image2]\n";
  cerr << "    -> Mesure the sum of squared differences between two images\n";
  cerr << "\n";
  cerr << "-Surf_Av_Std_values [RefIma] [Thresh] [Values]\n";
  cerr << "    -> Mesure the average and std of the intensities in [Values] at the surface of [RefIma] at the threshold [Thresh]\n";
  exit(1);
}

int main(int argc, char **argv){
  bool done;
  int int1,int2;
  char File1[256];
  char File2[256];
  char File3[256];
  char File4[256];
  char File5[256];
  
  done=false;
  
  // Check command line
  if (argc < 2) {
    usage();
  }
  else {
    //+++++++++++++++++++ image transport +++++++++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-transport") == 0) {
      argc--; argv++;
      int1 = atoi(argv[1]); //first time step
      argc--; argv++;
      strcpy(File1,argv[1]); //image to transport
      argc--; argv++;
      strcpy(File2,argv[1]); //velocity field X
      argc--; argv++;
      strcpy(File3,argv[1]); //velocity field Y
      argc--; argv++;
      strcpy(File4,argv[1]); //velocity field Z
      argc--; argv++;
      int2 = atoi(argv[1]); //final time step
      argc--; argv++;
      strcpy(File5,argv[1]); //transported image
      argc--; argv++;
      
      transportImage(File1,File5,int1,int2,File2,File3,File4);
      done = true;
    }
    //+++++++++++++   transport several images +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-MultiTransport") == 0) {
      MultiTransport(argc,argv);
      done = true;
    }
    //+++++++++++++   image weighted sum  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-weightedSum") == 0) {
      weightedSum(argc,argv);
      done = true;
    }
    //+++++++++++++   images weighted mean  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-weightedMean") == 0) {
      weightedMean(argc,argv);
      done = true;
    }
    //+++++++++++++   images weighted standard deviation  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-weightedStd") == 0) {
      weightedStd(argc,argv);
      done = true;
    }
    //+++++++++++++   sum of squared differences between two images  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-SSD") == 0) {
      MeasureSSD(argc,argv);
      done = true;
    }
    //+++++++++++++   surfacic average and std values  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Surf_Av_Std_values") == 0) {
      Surf_Av_Std_values(argc,argv);
      done = true;
    }

  }
  return 0;
}
