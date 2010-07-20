/*=========================================================================
  Date      : $Date: 12.02.2010$
=========================================================================*/

#include <irtkLargeDeformationSciCalcPack.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           1:   FUNCTIONS FOR THE CLASS "ScalarField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///constructor
ScalarField::ScalarField(void){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  this->NT=0;
}

///destructor
ScalarField::~ScalarField(void){}


///put a value
void ScalarField::P(float value,int x,int y,int z,int t){
  this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x]=value;
}

/// add a value
void ScalarField::Add(float value,int x,int y,int z,int t){
  this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x]+=value;
}

///put a the same value at every points of the scalar field
void ScalarField::PutToAllVoxels(float cste,int t)
{
  int x,y,z;
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) { this->P(cste,x,y,z,t); }
}

///get a value
float ScalarField::G(int x,int y,int z,int t){
  return this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x];
}

///get a value using linear interpolation
float ScalarField::G(float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  //values out of the image
  if (x<0.) x=0.0001;
  if (x>=this->NX-1.) x=this->NX-1.0001;
  if (y<0.) y=0.0001;
  if (y>=this->NY-1.) y=this->NY-1.0001;
  if (z<0.) z=0.0001;
  if (z>=this->NZ-1.) z=this->NZ-1.0001;
  if (t<0) t=0;
  if (t>this->NT-1) t=this->NT-1;
  
  //closest entire value
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  //interpolation
  if (this->NZ==1){ //2D IMAGE
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
    
    InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
  }
  else{//3D IMAGE
    wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
    wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
    
    InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
  }
  
  return InterpoGreyLevel;
}

///same as above
float ScalarField::G(double x,double y,double z,int t){
  return this->G((float)x,(float)y,(float)z,t);
}

///get the maximum absolute values out of the scalar field
float ScalarField::GetMaxAbsVal(int t){
  float max=0.0;
  int x,y,z;
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
  {
    if(max<abs(this->G(x,y,z,t))){max = abs(this->G(x,y,z,t));}
  }
  return max;
}


///read a scalar field (in a nifti image) -> DEPENDANCE TO IRTK
void ScalarField::Read(char * ImageName){
  irtkGenericImage<float> InputImage;
  int x,y,z,t;
  
  //read the image in the irtk format
  InputImage.Read(ImageName);
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=InputImage.GetX())||(this->NY!=InputImage.GetY())||(this->NZ!=InputImage.GetZ())||(this->NT!=InputImage.GetT()))
      cout << "WARNING: THE SIZE OF A NON-NULL SCALAR FIELD IS CHANGED\n";
  
  
  //fill the parameters of the class and allocate the memory for the image
  this->NX=InputImage.GetX();
  this->NY=InputImage.GetY();
  this->NZ=InputImage.GetZ();
  this->NT=InputImage.GetT();
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  //cast the image to the format used in the class
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          this->P(static_cast<float>(InputImage.Get(x, y, z, t)),x,y,z,t);
}

///read a scalar field and perform linear interpolation to give it a specific size
void ScalarField::Read_and_Interpolate(char * ImageName,int NBX,int NBY,int NBZ){
  ScalarField OrigSF;
  int x,y,z,t;
  float x2,y2,z2;
  
  //read the scalar field at the original format
  OrigSF.Read(ImageName);
  
  //fill the parameters of the class and allocate the memory for the scalar field
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=OrigSF.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  //interpolate the original image
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t),x,y,z,t);
        }
      }
    }
  }
}


///create a void scalar field. All the values are initialize to 'cste' which is null by default
void ScalarField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT,float cste){
  int x,y,z,t;
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=NBX)||(this->NY!=NBY)||(this->NZ!=NBZ)||(this->NT!=NBT))
      cout << "WARNING: THE SIZE OF A NON-NULL SCALAR FIELD IS CHANGED\n";
  
  //image size
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=NBT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  
  //allocate memory to cast (and eventually transform) the original template and target images
  //    -->  ScalarField[ptSF(x,y,z)]= gray level at (x,y,z)
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  //set all entries of the field at 0.
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          this->P(cste,x,y,z,t);
}


///write a scalar field (from a nifti image) -> DEPENDENCE TO IRTK
void ScalarField::Write(char * OutputImageName, char * ImageForHeaderName){
  irtkRealImage OutputImage;
  irtkRealImage ImageForHeader;
  irtkImageAttributes ImageAttribs;
  int x,y,z,t;
  
  //read the image from which we will extract the header
  ImageForHeader.Read(ImageForHeaderName);
  
  //extract and transform the header
  ImageAttribs=ImageForHeader.GetImageAttributes();
  ImageAttribs._x=this->NX;
  ImageAttribs._y=this->NY;
  ImageAttribs._z=this->NZ;
  ImageAttribs._t=this->NT;
  
  
  //cast the image form the format used in the class to the one of irtk
  OutputImage = irtkGenericImage<float>(ImageAttribs);

  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          OutputImage.Put(x, y, z, t,this->G(x,y,z,t));
  
  //write the image
  OutputImage.Write(OutputImageName);
}
  
///return the number of voxels in a ScalarField
int ScalarField::GetNbVoxels(){
  return this->NX *this->NY*this->NZ*this->NT;
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           2:   FUNCTIONS FOR THE CLASS "VectorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///constructor
VectorField::VectorField(void){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  this->NT=0;
}

///destructor
VectorField::~VectorField(void){}

///put a value
void VectorField::P(float value,int IdDirec,int x,int y,int z,int t){
  this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x]=value;
}
/// add a value
void VectorField::Add(float value,int IdDirec,int x,int y,int z,int t){
  this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x]+=value;
}

///put the same value at all entries of the vector field
void VectorField::PutToAllVoxels(float cste,int t){
  int i,x,y,z;
  for (i=0;i<3;i++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) { this->P(cste,i,x,y,z,t); }
}


///get a value
float VectorField::G(int IdDirec,int x,int y,int z,int t){
  return this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x];
}

///get a value using linear interpolation
float VectorField::G(int IdDirec,float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  
  //values out of the image
  if (x<0.) x=0.0001;
  if (x>=this->NX-1.) x=this->NX-1.0001;
  if (y<0.) y=0.0001;
  if (y>=this->NY-1.) y=this->NY-1.0001;
  if (z<0.) z=0.0001;
  if (z>=this->NZ-1.) z=this->NZ-1.0001;
  if (t<0) t=0;
  if (t>this->NT-1) t=this->NT-1;
  
  //closest entire value
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  //interpolation
  if (this->NZ==1){ //2D IMAGE
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
    
    InterpoGreyLevel= wmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
  }
  else{//3D IMAGE
    wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
    wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
    
    InterpoGreyLevel= wmmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wmpp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
    InterpoGreyLevel+=wpmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wpmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    InterpoGreyLevel+=wppp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
  }
  
  return InterpoGreyLevel;
}

///same as above
float VectorField::G(int IdDirec,double x,double y,double z,int t){
  return this->G(IdDirec,(float) x,(float) y,(float) z,t);
}

///get the maximum of the absolute values of the vector field
float VectorField::GetMaxAbsVal(int t)
{
  float max=0.0;
  int direc,x,y,z;
  for(direc=0;direc<3;direc++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
  {
    if(max<abs(this->G(direc,x,y,z,t))){max = abs(this->G(direc,x,y,z,t));}
  }
  return max;
}


///read a vector field (in 3 nifti images -> X, Y, Z) -> DEPENDENCE TO IRTK
void VectorField::Read(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z){
  irtkGenericImage<float> VecField_X;
  irtkGenericImage<float> VecField_Y;
  irtkGenericImage<float> VecField_Z;
  int x,y,z,t;
  
  //read the vector fields in the irtk format
  VecField_X.Read(NameVecField_X);
  VecField_Y.Read(NameVecField_Y);
  VecField_Z.Read(NameVecField_Z);
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=VecField_X.GetX())||(this->NY!=VecField_X.GetY())||(this->NZ!=VecField_X.GetZ())||(this->NT!=VecField_X.GetT()))
      cout << "WARNING: THE SIZE OF A NON-NULL VECTOR FIELD IS CHANGED\n";
  
  //fill the parameters of the class and allocate the memory for the image
  this->NX=VecField_X.GetX();
  this->NY=VecField_X.GetY();
  this->NZ=VecField_X.GetZ();
  this->NT=VecField_X.GetT();
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXtYtZtT=this->NXtYtZ*this->NT;
  this->VecField = new float [this->NXtYtZtT*3];
  
  //cast the image to the format used in the class
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          this->P(static_cast<float>(VecField_X.Get(x, y, z, t)),0,x,y,z,t);
  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          this->P(static_cast<float>(VecField_Y.Get(x, y, z, t)),1,x,y,z,t);
  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          this->P(static_cast<float>(VecField_Z.Get(x, y, z, t)),2,x,y,z,t);
}


///read a scalar vector and perform linear interpolation to give it a specific size
void VectorField::Read_and_Interpolate(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z,int NBX,int NBY,int NBZ,int rescaleVF){
  ScalarField OrigSF;
  int x,y,z,t;
  float x2,y2,z2;
  float scaleFactor;
  
  //0) init
  scaleFactor=1.;
  
  //1) X DIRECTION
  //1.1) read the scalar field in direction X at the original format
  OrigSF.Read(NameVecField_X);
  
  //1.2) fill the parameters of the class and allocate the memory for the vector field
  //(the directions Y and Z are supposed in the same format)
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=OrigSF.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXtYtZtT=this->NXtYtZ*this->NT;
  this->VecField = new float [this->NXtYtZtT*3];
  
  //1.3) manage the scale factor
  if (rescaleVF!=0) scaleFactor=((float)this->NX)/((float)OrigSF.NX);
  
  //1.4) interpolate the original image to compute the vector field in direction X
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,0,x,y,z,t);
        }
      }
    }
  }
  
  //2) Y DIRECTION
  //2.1) read the scalar field in direction Y
  OrigSF.Read(NameVecField_Y);
  
  //2.2) manage the scale factor
  if (rescaleVF!=0) scaleFactor=((float)this->NY)/((float)OrigSF.NY);
  
  //2.3) interpolate the original image to compute the vector field in direction Y
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,1,x,y,z,t);
        }
      }
    }
  }
  
  
  //3) Z DIRECTION
  //3.1) read the scalar field in direction Z
  OrigSF.Read(NameVecField_Z);
  
  //3.2) manage the scale factor
  if (rescaleVF!=0) scaleFactor=((float)this->NZ)/((float)OrigSF.NZ);
  
  //3.3) interpolate the original image to compute the vector field in direction Z
  for(t=0;t<this->NT;t++){
    for(z=0;z<this->NZ;z++){ 
      z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
      for(y=0;y<this->NY;y++){ 
        y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
        for(x=0;x<this->NX;x++){
          x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
          this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,2,x,y,z,t);
        }
      }
    }
  }
}



///constructor
void VectorField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT){
  int x,y,z,t,direc;
  
    //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=NBX)||(this->NY!=NBY)||(this->NZ!=NBZ)||(this->NT!=NBT))
      cout << "WARNING: THE SIZE OF A NON-NULL VECTOR FIELD IS CHANGED\n";
  
  //image size
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=NBT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXtYtZtT=this->NXtYtZ*this->NT;
  
  //allocate memory to cast (and eventually transform) the original template and target images
  //    -->  ScalarField[ptSF(x,y,z)]= gray level at (x,y,z)
  this->VecField = new float [this->NXtYtZtT*3];

  //fill the image with 0.
  for(direc=0;direc<3;direc++) for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
            this->P(0.,direc,x,y,z,t);
}

///write a vector field (from 3 nifti images -> X, Y, Z) -> DEPENDENCE TO IRTK
void VectorField::Write(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z){
  irtkRealImage OutputImage;
  int x,y,z,t;

  //irtk image to cast the vector fields form the format used in the class to the one of irtk
  OutputImage = irtkGenericImage<float>(this->NX, this->NY, this->NZ,this->NT);
  
  //write the vector field in X direction
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          OutputImage.Put(x, y, z, t,this->G(0,x,y,z,t));
  
  OutputImage.Write(NameVecField_X);
  
  //write the vector field in Y direction
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          OutputImage.Put(x, y, z, t,this->G(1,x,y,z,t));
  
  OutputImage.Write(NameVecField_Y);
  
  //write the vector field in Z direction
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
          OutputImage.Put(x, y, z, t,this->G(2,x,y,z,t));
  
  OutputImage.Write(NameVecField_Z);
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           3: CLASS TO PERFORM CONVOLUTION AND DECONVOLUTION USING FFT
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///Constructor
FFTconvolver3D::FFTconvolver3D(){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  
  this->NXfft=0;
  this->NYfft=0;
  this->NZfft=0;
}

///destructor
FFTconvolver3D::~FFTconvolver3D(void){}


///Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
///4 Gaussians (set some weights to 0 if less Gaussians are required)
///* NX, NY, NZ: is the size of the input image
///* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
///* w2,sX2,sY2,sZ2,: weight of the 2nd Gaussian kernel and std. dev. in direction X, Y, Z
///* w3,sX3,sY3,sZ3,: weight of the 3rd Gaussian kernel and std. dev. in direction X, Y, Z
///* w4,sX4,sY4,sZ4,: weight of the 4th Gaussian kernel and std. dev. in direction X, Y, Z
void FFTconvolver3D::InitiateConvolver(int NBX,int NBY, int NBZ,float w1,float sX1,float sY1,float sZ1,float w2,float sX2,float sY2,float sZ2,float w3,float sX3,float sY3,float sZ3,float w4,float sX4,float sY4,float sZ4){
  
  //set the size of the original image
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  
  //set the size of images for the FFT
  this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
  this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
  this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
  
  cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << "\n";
  
  //allocate memory for the images for the FFT
  this->RealSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - real part
  this->ImagSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - imaginary part
  this->RealFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - real part
  this->ImagFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - imaginary part
  
  //allocate memory for the temporary image
  this->ImageTemp.CreateVoidField(this->NXfft,this->NYfft,this->NZfft);
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf4AnisotropicGaussianFilters(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4);
}

///design a kernel that is the sum of up to 4 Gaussians and transform it in Fourier spaces
void FFTconvolver3D::MakeSumOf4AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4){
  int k,x,y,z;
  float SumLoc;
  float weight,sigmaX,sigmaY,sigmaZ;
  
  //compute and save the 4 kernels
  for (k=0;k<4;k++){
    //parameters of the current kernel
    if (k==3)     {weight=weight4; sigmaX=sigmaX4; sigmaY=sigmaY4; sigmaZ=sigmaZ4;}
    else if (k==2){weight=weight3; sigmaX=sigmaX3; sigmaY=sigmaY3; sigmaZ=sigmaZ3;}
    else if (k==1){weight=weight2; sigmaX=sigmaX2; sigmaY=sigmaY2; sigmaZ=sigmaZ2;}
    else          {weight=weight1; sigmaX=sigmaX1; sigmaY=sigmaY1; sigmaZ=sigmaZ1;}
    
    //design the current kernel with no influence of the weight
    for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++)
          this->ImageTemp.P((float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
          this->ImageTemp.P((float)(exp( -(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++)
          this->ImageTemp.P((float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
          this->ImageTemp.P((float)(exp( -(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++)
          this->ImageTemp.P((float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
          this->ImageTemp.P((float)(exp(-(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++)
          this->ImageTemp.P((float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
          this->ImageTemp.P((float)(exp(-(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    //normalization of the current filter and copy in RealFilterForFFT
    SumLoc=0.;
    for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) SumLoc+=this->ImageTemp.G(x,y,z);
    if (k==0){
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(weight*this->ImageTemp.G(x,y,z)/SumLoc,x,y,z);
    }
    else{
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(this->RealFilterForFFT.G(x,y,z)+weight*this->ImageTemp.G(x,y,z)/SumLoc,x,y,z);
    }
  }
  
  //set ImagFilterForFFT to 0 in case it contains something
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT.P(0.,x,y,z);
  
  //Transform RealFilterForFFT and ImagFilterForFFT in Fourier spaces
  this->DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
}


///Fast Fourier Transform of numerical recipies (slighly modified)
void FFTconvolver3D::four1NR(float data[], unsigned long nn, int isign){
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2){
    if (j>i){
      tempr=data[j]; data[j]=data[i]; data[i]=tempr;
      tempr=data[j+1]; data[j+1]=data[i+1]; data[i+1]=tempr;
    }
    m=n >> 1;
    while ((m>=2) && (j>m)){
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

///change the kernel of the convolver (same notations as the constructor)
void FFTconvolver3D::ChangeKernel(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4){
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf4AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,weight2,sigmaX2,sigmaY2,sigmaZ2,weight3,sigmaX3,sigmaY3,sigmaZ3,weight4,sigmaX4,sigmaY4,sigmaZ4);
}

///Fast Fourier Transform
void FFTconvolver3D::DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));

  
  //2) perform the fft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->G(x, y, z);
      dataX[2*x+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataX, (unsigned long)SizeX, 1);
    for (x = 0; x < SizeX; x++){
      RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
      ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX, x, y, z);
    }
  }
  delete dataX;
  
  //3) perform the fft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->G(x, y, z);
      dataY[2*y+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataY, (unsigned long)SizeY, 1);
    for (y = 0; y < SizeY; y++){
      RealSignal->P(dataY[2*y+1]/SqrtSizeY,x, y, z);
      ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
    }
  }
  delete dataY;
  
  
  //4) perform the fft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->G(x, y, z);
      dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataZ, (unsigned long)SizeZ, 1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->P(dataZ[2*z+1]/SqrtSizeZ,x, y, z);
      ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ, x, y, z);
    }
  }
  delete dataZ;
}


///Inverse Fast Fourier Transform
void FFTconvolver3D::InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  
  //2) perform the ifft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->G(x, y, z);
      dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataZ, (unsigned long)SizeZ, -1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->P(dataZ[2*z+1]/SqrtSizeZ, x, y, z);
      ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ,x, y, z);
    }
  }
  delete dataZ;
  
  //3) perform the ifft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->G(x, y, z);
      dataY[2*y+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataY, (unsigned long)SizeY, -1);
    for (y = 0; y < SizeY; y++){
      RealSignal->P(dataY[2*y+1]/SqrtSizeY, x, y, z);
      ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
    }
  }
  delete dataY;
  
  //4) perform the ifft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->G(x, y, z);
      dataX[2*x+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataX, (unsigned long)SizeX, -1);
    for (x = 0; x < SizeX; x++){
      RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
      ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX,x, y, z);
    }
  }
  delete dataX;
}


///convolution of a 3D scalar field using the predifined kernel
void FFTconvolver3D::Convolution(ScalarField * SF){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  //1) Copy the orginal image in the image that will be transformed
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
  
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
    this->RealSignalForFFT.P(SF->G(x,y,z),x,y,z);
  }
  
  //2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
  this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //3) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
  
  for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
    a=this->RealSignalForFFT.G(x, y, z);
    b=this->ImagSignalForFFT.G(x, y, z);
    c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
    d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
    
    this->RealSignalForFFT.P(a*c-b*d, x, y, z);
    this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
  }
  
  //4) IFFT
  this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //5) Copy the image that has been convolved in the orginal image
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
    SF->P(this->RealSignalForFFT.G(x,y,z),x,y,z);
  }
}


///convolution of the real scalar field defined inside of the class
void FFTconvolver3D::Convolution(){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  //1) Set to 0. all values that cannot be accessed by outside of the class
  for (z=this->NZ;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++)        for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZ;z++)           for (y=this->NY;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZ;z++)           for (y=0;y<this->NY;y++)           for (x=this->NX;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) 
        this->ImagSignalForFFT.P(0.,x,y,z);
  
  //2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
  this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //3) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
  
  for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
    a=this->RealSignalForFFT.G(x, y, z);
    b=this->ImagSignalForFFT.G(x, y, z);
    c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
    d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
    
    this->RealSignalForFFT.P(a*c-b*d, x, y, z);
    this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
  }
  
  //4) IFFT
  this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
}

///put a value in the real part of the field that is transformed by the class
void FFTconvolver3D::P(float value,int x,int y, int z){
  this->RealSignalForFFT.P(value,x,y,z);
}

///put a value in the real part of the field that is transformed by the class
float FFTconvolver3D::G(int x,int y, int z){
  return this->RealSignalForFFT.G(x,y,z);
}

///deconvolution of a 3D scalar field using the predifined kernel
/// !!! NOT VALIDATED !!!
void FFTconvolver3D::Deconvolution(ScalarField * SF){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  cout << "DECONVOLUTION SHOULD BE USED CARREFULLY HERE - NOT VALIDATED!!!\n";
  
  //1) Copy the orginal image in the image that will be transformed
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
  
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
    this->RealSignalForFFT.P(SF->G(x,y,z),x,y,z);
  }
  //2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
  this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);

  //3) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
  
  for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
    a=this->RealSignalForFFT.G(x, y, z);
    b=this->ImagSignalForFFT.G(x, y, z);
    c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
    d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
    
    this->RealSignalForFFT.P((a*c+b*d)/(c*c+d*d), x, y, z);
    this->ImagSignalForFFT.P((c*b-a*d)/(c*c+d*d),x, y, z);
  }
  //4) IFFT
  this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
  
  //5) Copy the image that has been deconvolved in the orginal image
  for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
    SF->P(this->RealSignalForFFT.G(x,y,z),x,y,z);
  }
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           4: LOW LEVEL FUNCTIONS MAKING USE OF THE CLASSES ScalarField AND VectorField 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
void Cpt_Grad_ScalarField(ScalarField * SField,VectorField * Gradient,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  NBX=SField->NX;
  NBY=SField->NY;
  NBZ=SField->NZ;
  NBT=SField->NT;
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
  
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)){
      Gradient->CreateVoidField(NBX,NBY,NBZ);
      cout << "Gradient added in Cpt_Grad_ScalarField\n";
    }
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) gradient in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z);
      Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z);
      Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z);
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
      
      //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0);
      Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)||(NBT!=Gradient->NT))
      Gradient->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //gradient in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z,t);
        Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z,t);
        Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z,t);
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0,t);
        Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0,t);
      }
    }
  }
}

//Compute (d VField(X) / d x) + (d VField(Y) / d y) + (d VField(Z) / d z) and put the result in 'GradScalVF'
//where 'VField' is a vector field and 'GradScalVF' a scalar field
void Cpt_Grad_Scal_VectorField(VectorField * VField,ScalarField * GradScalVF,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float GradX,GradY,GradZ;
  
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  NBT=VField->NT;
  
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
  
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=GradScalVF->NX)||(NBY!=GradScalVF->NY)||(NBZ!=GradScalVF->NZ))
      GradScalVF->CreateVoidField(NBX,NBY,NBZ);
    
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) sum of gradients in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      GradX=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
      GradY=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
      GradZ=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
      GradScalVF->P(GradX+GradY+GradZ,x,y,z);
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      GradX=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
      GradY=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
      GradScalVF->P(GradX+GradY,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=GradScalVF->NX)||(NBY!=GradScalVF->NY)||(NBZ!=GradScalVF->NZ)||(NBT!=GradScalVF->NT))
      GradScalVF->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //sum of gradients in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        GradX=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
        GradY=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
        GradZ=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
        GradScalVF->P(GradX+GradY+GradZ,x,y,z,t);
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        GradX=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
        GradY=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
        GradScalVF->P(GradX+GradY,x,y,0,t);
      }
    }
  }
}


//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;

  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  NBT=VField->NT;
  

  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
  
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ))
      DetJ->CreateVoidField(NBX,NBY,NBZ);
    
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) sum of gradients in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
      d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
      d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
      d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
      d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
      d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
      DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z);
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
      DetJ->P(d11*d22-d21*d12,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ)||(NBT!=DetJ->NT))
      DetJ->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //sum of gradients in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
        d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
        d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
        d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
        d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
        d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
        d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
        d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
        d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
        DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z,t);
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
        d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
        d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
        d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
        DetJ->P(d11*d22-d21*d12,x,y,0,t);
      }
    }
  }
}

//Compute the 4D forward mapping 'Fmap' from the velocity field 'VeloField'
void ForwardMappingFromVelocityField(VectorField * VeloField,VectorField * Fmap,VectorField * MapSrcImag,int ConvergenceSteps,float DeltaX,int StartPoint){
  float VecTemp[3];
  float VecTemp2[3];
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  
  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
    //allocate memory in GradScalVF if not done
  if ((NBX!=Fmap->NX)||(NBY!=Fmap->NY)||(NBZ!=Fmap->NZ)||(NBT!=Fmap->NT))
    Fmap->CreateVoidField(NBX,NBY,NBZ,NBT);

  
  //JO at the first time subdivision
  if ((StartPoint<0)||(StartPoint>NBT-1)){ //normal ForwardMapping from the initial time
    StartPoint=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      Fmap->P(MapSrcImag->G(0,x,y,z),0,x,y,z,0);
      Fmap->P(MapSrcImag->G(1,x,y,z),1,x,y,z,0);
      Fmap->P(MapSrcImag->G(2,x,y,z),2,x,y,z,0);
    }
  }
  else{
    for (t=0;t<=StartPoint;t++){  //ForwardMapping from the a time between 0 and 1
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        Fmap->P(MapSrcImag->G(0,x,y,z),0,x,y,z,t);
        Fmap->P(MapSrcImag->G(1,x,y,z),1,x,y,z,t);
        Fmap->P(MapSrcImag->G(2,x,y,z),2,x,y,z,t);
      }
    }
  }
  
  
  //JO at the other time subdivisions
  for (t=StartPoint+1;t<NBT;t++){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        //find the original coordinates
        Fmap->P(Fmap->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Fmap->P(Fmap->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Fmap->P(Fmap->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[1]=VeloField->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[2]=VeloField->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Fmap->P(Fmap->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Fmap->P(Fmap->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Fmap->P(Fmap->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
  }
}


//Compute the 4D backward mapping 'Bmap' from the velocity field 'VeloField'
void BackwardMappingFromVelocityField(VectorField * VeloField,VectorField * Bmap,VectorField * MapTrgImag,int ConvergenceSteps,float DeltaX,int StartPoint){
  float VecTemp[3];
  float VecTemp2[3];
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;

  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //allocate memory in GradScalVF if not done
  if ((NBX!=Bmap->NX)||(NBY!=Bmap->NY)||(NBZ!=Bmap->NZ)||(NBT!=Bmap->NT))
    Bmap->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //Bmap at the last time subdivision
  if ((StartPoint<0)||(StartPoint>NBT-1)){ //normal ForwardMapping from the initial time
    StartPoint=NBT-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      Bmap->P(MapTrgImag->G(0,x,y,z),0,x,y,z,NBT-1);
      Bmap->P(MapTrgImag->G(1,x,y,z),1,x,y,z,NBT-1);
      Bmap->P(MapTrgImag->G(2,x,y,z),2,x,y,z,NBT-1);
    }
  }
  else{
    for (t=NBT-1;t>=StartPoint;t--){  //ForwardMapping from the a time between 0 and 1
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        Bmap->P(MapTrgImag->G(0,x,y,z),0,x,y,z,t);
        Bmap->P(MapTrgImag->G(1,x,y,z),1,x,y,z,t);
        Bmap->P(MapTrgImag->G(2,x,y,z),2,x,y,z,t);
      }
    }
  }
  
  //Bmap at the other time subdivisions
  for (t=StartPoint-1;t>=0;t--){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        Bmap->P(Bmap->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Bmap->P(Bmap->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Bmap->P(Bmap->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[1]=VeloField->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[2]=VeloField->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Bmap->P(Bmap->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Bmap->P(Bmap->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Bmap->P(Bmap->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
  }
}


//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialFmap' which is the partial forward mapping due to the contribution of PartialVeloField.
//Imporant: Should only be used for the transportation of the source image
void ComputePartialForwardMapping(VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialFmap,VectorField * MapSrcImag,VectorField * MapTrgImag,int ConvergenceSteps,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  float x1,y1,z1;
  float x2,y2,z2;
  float x3,y3,z3;
  float x4,y4,z4;
  float xS,yS,zS;
  int t;
  float DeltaT_div_DeltaX;
  VectorField TotalBmap;  //total Backward mapping from the [new time_sub-1] to 0
  VectorField TotalBmap2;  //total Backward mapping from the [new time_sub] to 0

  //1) initialisation
  //1.1) constants
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in GradScalVF if not done
  if ((NBX!=PartialFmap->NX)||(NBY!=PartialFmap->NY)||(NBZ!=PartialFmap->NZ)||(NBT!=PartialFmap->NT))
    PartialFmap->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //1.3) allocate memory for TotalBmap
  TotalBmap.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalBmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //1.4) PartialFmap at the first time subdivision
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    PartialFmap->P(MapSrcImag->G(0,x,y,z),0,x,y,z,0);
    PartialFmap->P(MapSrcImag->G(1,x,y,z),1,x,y,z,0);
    PartialFmap->P(MapSrcImag->G(2,x,y,z),2,x,y,z,0);
  }
  
  //2) PartialFmap at the other time subdivisions
  for (t=1;t<NBT;t++){
    
    //2.1) compute the total backward mapping from t-1 to 0
    BackwardMappingFromVelocityField(VeloField,&TotalBmap,MapTrgImag,ConvergenceSteps,1,t-1);
    
    //2.2) compute the total backward mapping from t to 0
    BackwardMappingFromVelocityField(VeloField,&TotalBmap2,MapTrgImag,ConvergenceSteps,1,t);
    
    //2.3) compute partial forward map at t from the one at t-1
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //2.3.1) first estimation
      //2.3.1.a) first guess of where the information comes from at t-1
      x1=PartialFmap->G(0,x,y,z,t-1); y1=PartialFmap->G(1,x,y,z,t-1); z1=PartialFmap->G(2,x,y,z,t-1);
      x2=TotalBmap.G(0,x1,y1,z1,0);   y2=TotalBmap.G(1,x1,y1,z1,0);   z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
      
      xS=x-PartialVeloField->G(0,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      yS=y-PartialVeloField->G(1,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      zS=z-PartialVeloField->G(2,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
      
      //2.3.1.b) first transport of the information
      PartialFmap->P(PartialFmap->G(0,xS,yS,zS,t-1),0,x,y,z,t);
      PartialFmap->P(PartialFmap->G(1,xS,yS,zS,t-1),1,x,y,z,t);
      PartialFmap->P(PartialFmap->G(2,xS,yS,zS,t-1),2,x,y,z,t);
      
      
      //2.3.1) leap frog style improvement of the estimation
      for (i=0;i<ConvergenceSteps*2;i++){
        //2.3.2.a) where the information comes from at t-1
        x1=PartialFmap->G(0,xS,yS,zS,t-1); y1=PartialFmap->G(1,xS,yS,zS,t-1); z1=PartialFmap->G(2,xS,yS,zS,t-1);
        x2=TotalBmap.G(0,x1,y1,z1,0);      y2=TotalBmap.G(1,x1,y1,z1,0);      z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
        
        x3=PartialFmap->G(0,x,y,z,t);   y3=PartialFmap->G(1,x,y,z,t);   z3=PartialFmap->G(2,x,y,z,t);
        x4=TotalBmap2.G(0,x3,y3,z3,0);  y4=TotalBmap2.G(1,x3,y3,z3,0);  z4=TotalBmap2.G(2,x3,y3,z3,0); //TotalBmap2 -> t
        
        xS=x-(PartialVeloField->G(0,x2,y2,z2,t-1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        yS=y-(PartialVeloField->G(1,x2,y2,z2,t-1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        zS=z-(PartialVeloField->G(2,x2,y2,z2,t-1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        
        //2.3.1.b) update the transport of the information
        PartialFmap->P(PartialFmap->G(0,xS,yS,zS,t-1),0,x,y,z,t);
        PartialFmap->P(PartialFmap->G(1,xS,yS,zS,t-1),1,x,y,z,t);
        PartialFmap->P(PartialFmap->G(2,xS,yS,zS,t-1),2,x,y,z,t);
      }
    }
  }
}


//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField). MapTrgImag is the mapping of Trg image (possibly not the identity).
//This function then computes 'PartialBmap' which is the partial backward mapping due to the contribution of PartialVeloField.
//Imporant: Should only be used for the transportation of the target image
void ComputePartialBackwardMapping(VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialBmap,VectorField * MapSrcImag,int ConvergenceSteps,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  float x1,y1,z1;
  float x2,y2,z2;
  float x3,y3,z3;
  float x4,y4,z4;
  float xS,yS,zS;
  int t;
  float DeltaT_div_DeltaX;
  VectorField TotalFmap;  //total forward mapping from the [new time_sub-1] to 0
  VectorField TotalFmap2;  //total forwrd mapping from the [new time_sub] to 0

  //1) initialisation
  //1.1) constants
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in GradScalVF if not done
  if ((NBX!=PartialBmap->NX)||(NBY!=PartialBmap->NY)||(NBZ!=PartialBmap->NZ)||(NBT!=PartialBmap->NT))
    PartialBmap->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //1.3) allocate memory for TotalFmap
  TotalFmap.CreateVoidField(NBX,NBY,NBZ,NBT);
  TotalFmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
  
  //1.4) PartialBmap at the first time subdivision
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    PartialBmap->P(MapSrcImag->G(0,x,y,z),0,x,y,z,NBT-1);
    PartialBmap->P(MapSrcImag->G(1,x,y,z),1,x,y,z,NBT-1);
    PartialBmap->P(MapSrcImag->G(2,x,y,z),2,x,y,z,NBT-1);
  }
  
  //2) PartialBmap at the other time subdivisions
  for (t=NBT-2;t>=0;t--){
    
    //2.1) compute the total forward mapping from 0 to t+1
    ForwardMappingFromVelocityField(VeloField,&TotalFmap,MapSrcImag,ConvergenceSteps,1,t+1);
    
    //2.2) compute the total forward mapping from 0 to t
    ForwardMappingFromVelocityField(VeloField,&TotalFmap2,MapSrcImag,ConvergenceSteps,1,t);
    
    //2.3) compute partial forward map at t from the one at t+1
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //2.3.1) first estimation
      //2.3.1.a) first guess of where the information comes from at t-1
      //if ((y<1)||(y>38)) cout << "in3\n";

      x1=PartialBmap->G(0,x,y,z,t+1); y1=PartialBmap->G(1,x,y,z,t+1); z1=PartialBmap->G(2,x,y,z,t+1);
      x2=TotalFmap.G(0,x1,y1,z1,NBT-1);   y2=TotalFmap.G(1,x1,y1,z1,NBT-1);   z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t+1
      //if ((y<1)||(y>38)) cout << "in4\n";
      
      xS=x+PartialVeloField->G(0,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      yS=y+PartialVeloField->G(1,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      zS=z+PartialVeloField->G(2,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
      //if ((y<1)||(y>38)) cout << "in5\n";
      
      //2.3.1.b) first transport of the information
      PartialBmap->P(PartialBmap->G(0,xS,yS,zS,t+1),0,x,y,z,t);
      PartialBmap->P(PartialBmap->G(1,xS,yS,zS,t+1),1,x,y,z,t);
      PartialBmap->P(PartialBmap->G(2,xS,yS,zS,t+1),2,x,y,z,t);
      
      //2.3.1) leap frog style improvement of the estimation
      for (i=0;i<ConvergenceSteps*2;i++){
        //2.3.2.a) where the information comes from at t-1
        x1=PartialBmap->G(0,xS,yS,zS,t+1); y1=PartialBmap->G(1,xS,yS,zS,t+1); z1=PartialBmap->G(2,xS,yS,zS,t+1);
        x2=TotalFmap.G(0,x1,y1,z1,NBT-1);      y2=TotalFmap.G(1,x1,y1,z1,NBT-1);      z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t-1
        
        x3=PartialBmap->G(0,x,y,z,t);   y3=PartialBmap->G(1,x,y,z,t);   z3=PartialBmap->G(2,x,y,z,t);
        x4=TotalFmap2.G(0,x3,y3,z3,NBT-1);  y4=TotalFmap2.G(1,x3,y3,z3,NBT-1);  z4=TotalFmap2.G(2,x3,y3,z3,NBT-1); //TotalFmap2 -> t
        
        xS=x+(PartialVeloField->G(0,x2,y2,z2,t+1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        yS=y+(PartialVeloField->G(1,x2,y2,z2,t+1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        zS=z+(PartialVeloField->G(2,x2,y2,z2,t+1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
        
        //2.3.1.b) update the transport of the information
        PartialBmap->P(PartialBmap->G(0,xS,yS,zS,t+1),0,x,y,z,t);
        PartialBmap->P(PartialBmap->G(1,xS,yS,zS,t+1),1,x,y,z,t);
        PartialBmap->P(PartialBmap->G(2,xS,yS,zS,t+1),2,x,y,z,t);
      }
    }
  }
}



//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialLocBmap' which is the partial forward mapping ONLY AT 'TargetSubdiv' FROM 'SourceSubdiv' due to the contribution of PartialVeloField.
//-> PartialLocBmap represents where will be the coordinates of the points of time subdivision 'SourceSubdiv' after transportation to time subdivision 'TargetSubdiv'
//Important: no mapping of the source or target map is supported here
//Imporant: Should only be used for the transportation of the source image
void ComputeLagrangianPartialBackwardMapping(VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialLocBmap,int SourceSubdiv,int TargetSubdiv,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float x1,y1,z1;
  float x2,y2,z2;
  float x3,y3,z3;
  float DeltaT_div_DeltaX;
  float DX,DY,DZ;
  float DX2,DY2,DZ2;
  
  //1) initialisation
  //1.1) constants
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in GradScalVF if not done
  if ((NBX!=PartialLocBmap->NX)||(NBY!=PartialLocBmap->NY)||(NBZ!=PartialLocBmap->NZ)||(NBT!=1))
    PartialLocBmap->CreateVoidField(NBX,NBY,NBZ,1);
  
  //2) Compute the transportation of the points of time SourceSubdiv to TargetSubdiv
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    
    //initial coordinates
    x1=x*1.;  y1=y*1.;  z1=z*1.;  //for the transportation in the complete velocity field
    x2=x*1.;  y2=y*1.;  z2=z*1.;  //for the transportation in the partial velocity field
    
    //transportation
    for (t=SourceSubdiv;t<TargetSubdiv;t++){
      x3=x1; y3=y1; z3=z1;
      
      DX=VeloField->G(0,x3,y3,z3,t)*DeltaT_div_DeltaX;
      DY=VeloField->G(1,x3,y3,z3,t)*DeltaT_div_DeltaX;
      DZ=VeloField->G(2,x3,y3,z3,t)*DeltaT_div_DeltaX;
      
      DX2=(VeloField->G(0,x3,y3,z3,t)+VeloField->G(0,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
      DY2=(VeloField->G(1,x3,y3,z3,t)+VeloField->G(1,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
      DZ2=(VeloField->G(2,x3,y3,z3,t)+VeloField->G(2,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
      DX=DX2;
      DY=DY2;
      DZ=DZ2;
      
      x1=x1+DX;
      y1=y1+DY;
      z1=z1+DZ;
      
      x2=x2+(PartialVeloField->G(0,x3,y3,z3,t)+PartialVeloField->G(0,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
      y2=y2+(PartialVeloField->G(1,x3,y3,z3,t)+PartialVeloField->G(1,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
      z2=z2+(PartialVeloField->G(2,x3,y3,z3,t)+PartialVeloField->G(2,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
      
    }
    
    //save where will be the point x,y,z
    PartialLocBmap->P(x2,0,x,y,z);
    PartialLocBmap->P(y2,1,x,y,z);
    PartialLocBmap->P(z2,2,x,y,z);
  }
}


//Compute the projection of a 3D image 'ImageTime0' using Forward Mapping 'Fmap'.
//The image is projected at the time step 'TimeStepProj' of 'Fmap' and stored in
//'ImageTimeT'.
void ForwardProjectionOf3Dimage(ScalarField * ImageTime0,VectorField * Fmap,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z;
  
  for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    ImageTimeT->P(ImageTime0->G(Fmap->G(0,x,y,z,TimeStepProj), Fmap->G(1,x,y,z,TimeStepProj),Fmap->G(2,x,y,z,TimeStepProj)),x,y,z);
  }
}


//Compute the projection of a 3D image 'ImageTime0' using Backward Mapping 'Bmap'.
//The image is projected at the time step 'TimeStepProj' of 'Fmap' and stored in
//'ImageTimeT'.
void BackwardProjectionOf3Dimage(ScalarField * ImageTime1,VectorField * Bmap,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z;
  
  for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    ImageTimeT->P(ImageTime1->G(Bmap->G(0,x,y,z,TimeStepProj), Bmap->G(1,x,y,z,TimeStepProj),Bmap->G(2,x,y,z,TimeStepProj)),x,y,z);
  }
}



///By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
///'VeloField4Measure' in the length of the flow from each point of the field. The length of flow
///is returned in the 3D scalar field 'LengthOfFlow'
/// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
/// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
///   is computed.
void CptLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps,float DeltaX){
  float VecTemp[3];
  float VecTemp2[3];
  float NormVecTemp;
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  ScalarField PrevLengthOfFlow;
      
  //1) INITIALISATION
  //1.1) field size
  NBX=VeloField4Flow->NX;
  NBY=VeloField4Flow->NY;
  NBZ=VeloField4Flow->NZ;
  NBT=VeloField4Flow->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //1.2) allocate memory in LengthOfFlow at time t if not done
  if ((NBX!=LengthOfFlow->NX)||(NBY!=LengthOfFlow->NY)||(NBZ!=LengthOfFlow->NZ))
    LengthOfFlow->CreateVoidField(NBX,NBY,NBZ);
  
  //1.3) allocate memory of the PrevLengthOfFlow at time t+1
  PrevLengthOfFlow.CreateVoidField(NBX,NBY,NBZ);
  
  //1.4) PrevLengthOfFlow at the last time subdivision
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) PrevLengthOfFlow.P(0.,x,y,z);
  
  //JO at the other time subdivisions
  for (t=NBT-2;t>=0;t--){
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
      //init
      VecTemp[0]=0.; 
      VecTemp[1]=0.;
      VecTemp[2]=0.;
      
      //convergence
      for (i=0;i<ConvergenceSteps;i++){
        VecTemp2[0]=VeloField4Flow->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
        VecTemp2[1]=VeloField4Flow->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
        VecTemp2[2]=VeloField4Flow->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
        
        VecTemp[0]=(VecTemp2[0]+VeloField4Flow->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VecTemp2[1]+VeloField4Flow->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VecTemp2[2]+VeloField4Flow->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
      }
      
      //compute the lenght
      NormVecTemp =(float)pow((double)VeloField4Measure->G(0,x,y,z,t)*DeltaT_div_DeltaX,2.);
      NormVecTemp+=(float)pow((double)VeloField4Measure->G(1,x,y,z,t)*DeltaT_div_DeltaX,2.);
      NormVecTemp+=(float)pow((double)VeloField4Measure->G(2,x,y,z,t)*DeltaT_div_DeltaX,2.);
      NormVecTemp=sqrt(NormVecTemp);
      
      LengthOfFlow->P(NormVecTemp+PrevLengthOfFlow.G(x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),x,y,z);
    }
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
          PrevLengthOfFlow.P(LengthOfFlow->G(x,y,z),x,y,z);
  }
}



///compute the L_2 norm of the difference between two scalar fields
float CalcSqrtSumOfSquaredDif(ScalarField * I1,ScalarField * I2){
  int x,y,z;
  float L2_norm,tmp;
  
  L2_norm=0.;
  
  for (z=0;z<I1->NZ;z++) for (y=0;y<I1->NY;y++) for (x=0;x<I1->NX;x++){
    tmp=(I1->G(x,y,z)-I2->G(x,y,z));
    L2_norm+=tmp*tmp;
  }
  
  L2_norm=sqrt(L2_norm);
  
  return L2_norm;
}



///...
void TransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,int t)
{	int x,y,z;
	int NBX,NBY,NBZ;
        float d11,d12,d13,d21,d22,d23,d31,d32,d33;
	float temp;
	NBX=TempInvDiffeo->NX;
	NBY=TempInvDiffeo->NY;
	NBZ=TempInvDiffeo->NZ;
	for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
	  d11=(TempInvDiffeo->G(0,x+1,y,z,t)-TempInvDiffeo->G(0,x-1,y,z,t))/(2.*DeltaX);
	  d12=(TempInvDiffeo->G(0,x,y+1,z,t)-TempInvDiffeo->G(0,x,y-1,z,t))/(2.*DeltaX);
	  d13=(TempInvDiffeo->G(0,x,y,z+1,t)-TempInvDiffeo->G(0,x,y,z-1,t))/(2.*DeltaX);
	  d21=(TempInvDiffeo->G(1,x+1,y,z,t)-TempInvDiffeo->G(1,x-1,y,z,t))/(2.*DeltaX);
	  d22=(TempInvDiffeo->G(1,x,y+1,z,t)-TempInvDiffeo->G(1,x,y-1,z,t))/(2.*DeltaX);
	  d23=(TempInvDiffeo->G(1,x,y,z+1,t)-TempInvDiffeo->G(1,x,y,z-1,t))/(2.*DeltaX);
	  d31=(TempInvDiffeo->G(2,x+1,y,z,t)-TempInvDiffeo->G(2,x-1,y,z,t))/(2.*DeltaX);
	  d32=(TempInvDiffeo->G(2,x,y+1,z,t)-TempInvDiffeo->G(2,x,y-1,z,t))/(2.*DeltaX);
	  d33=(TempInvDiffeo->G(2,x,y,z+1,t)-TempInvDiffeo->G(2,x,y,z-1,t))/(2.*DeltaX);
	  temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0);
	  Momentum->P( temp* (d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13)),x,y,z);
	}
	//1.2.2) boundaries at 0.
	z=0;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	z=NBZ-1;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	y=0;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	y=NBY-1;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	x=0;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	x=NBX-1;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);

	//1.2.3) 2D image case
	//float max=0.0;
	if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
	  d11=(TempInvDiffeo->G(0,x+1,y,0,t)-TempInvDiffeo->G(0,x-1,y,0,t))/(2.*DeltaX);
	  d12=(TempInvDiffeo->G(0,x,y+1,0,t)-TempInvDiffeo->G(0,x,y-1,0,t))/(2.*DeltaX);
	  d21=(TempInvDiffeo->G(1,x+1,y,0,t)-TempInvDiffeo->G(1,x-1,y,0,t))/(2.*DeltaX);
	  d22=(TempInvDiffeo->G(1,x,y+1,0,t)-TempInvDiffeo->G(1,x,y-1,0,t))/(2.*DeltaX);
	  temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,0,t),TempInvDiffeo->G(1,x,y,0,t),TempInvDiffeo->G(2,x,y,0,t),0);
	  //if (max<abs(temp*(d11*d22-d21*d12))){max=abs(temp*(d11*d22-d21*d12));}
	  Momentum->P(temp*(d11*d22-d21*d12),x,y,0);
	}
	//cout << "c'est penible  "<< Momentum->GetMaxAbsVal() <<"\n";
}

/// Bug - Not VaLidated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void AddTransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,float cste, int t)
{	int x,y,z;
	int NBX,NBY,NBZ;
    float d11,d12,d13,d21,d22,d23,d31,d32,d33;
	float temp;
	NBX=TempInvDiffeo->NX;
	NBY=TempInvDiffeo->NY;
	NBZ=TempInvDiffeo->NZ;
	for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
	  d11=(TempInvDiffeo->G(0,x+1,y,z,t)-TempInvDiffeo->G(0,x-1,y,z,t))/(2.*DeltaX);
	  d12=(TempInvDiffeo->G(0,x,y+1,z,t)-TempInvDiffeo->G(0,x,y-1,z,t))/(2.*DeltaX);
	  d13=(TempInvDiffeo->G(0,x,y,z+1,t)-TempInvDiffeo->G(0,x,y,z-1,t))/(2.*DeltaX);
	  d21=(TempInvDiffeo->G(1,x+1,y,z,t)-TempInvDiffeo->G(1,x-1,y,z,t))/(2.*DeltaX);
	  d22=(TempInvDiffeo->G(1,x,y+1,z,t)-TempInvDiffeo->G(1,x,y-1,z,t))/(2.*DeltaX);
	  d23=(TempInvDiffeo->G(1,x,y,z+1,t)-TempInvDiffeo->G(1,x,y,z-1,t))/(2.*DeltaX);
	  d31=(TempInvDiffeo->G(2,x+1,y,z,t)-TempInvDiffeo->G(2,x-1,y,z,t))/(2.*DeltaX);
	  d32=(TempInvDiffeo->G(2,x,y+1,z,t)-TempInvDiffeo->G(2,x,y-1,z,t))/(2.*DeltaX);
	  d33=(TempInvDiffeo->G(2,x,y,z+1,t)-TempInvDiffeo->G(2,x,y,z-1,t))/(2.*DeltaX);
	  temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0);
	  Momentum->Add(cste* temp* (d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13)),x,y,z);
	}
	//1.2.2) boundaries at 0.
	z=0;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	z=NBZ-1;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	y=0;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	y=NBY-1;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	x=0;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	x=NBX-1;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);

	//1.2.3) 2D image case
	if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
	  d11=(TempInvDiffeo->G(0,x+1,y,0,t)-TempInvDiffeo->G(0,x-1,y,0,t))/(2.*DeltaX);
	  d12=(TempInvDiffeo->G(0,x,y+1,0,t)-TempInvDiffeo->G(0,x,y-1,0,t))/(2.*DeltaX);
	  d21=(TempInvDiffeo->G(1,x+1,y,0,t)-TempInvDiffeo->G(1,x-1,y,0,t))/(2.*DeltaX);
	  d22=(TempInvDiffeo->G(1,x,y+1,0,t)-TempInvDiffeo->G(1,x,y-1,0,t))/(2.*DeltaX);
	  temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,0,t),TempInvDiffeo->G(1,x,y,0,t));
	  Momentum->Add(cste*temp*(d11*d22-d21*d12),x,y,0);
	}
}

///...
void TransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image, int t)
{
	int x,y,z;
    for (z = 0; z < InitialImage->NZ; z++) for (y = 0; y < InitialImage->NY; y++) for (x = 0; x < InitialImage->NX; x++)
	{
		Image->P(InitialImage->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0),x,y,z,0);
	}
}

///...
void AddTransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image, float cste, int t)
{
	int x,y,z;
    for (z = 0; z < InitialImage->NZ; z++) for (y = 0; y < InitialImage->NY; y++) for (x = 0; x < InitialImage->NX; x++)
	{
		Image->Add(cste * InitialImage->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0),x,y,z,0);
	}
}

///...
void DeepCopy(VectorField *VectorField1,VectorField *VectorField2,int t)
{
	int i,x,y,z;
	for (i=0;i<3;i++)
	{
		for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
		{
			VectorField2->P(VectorField1->G(i,x,y,z),i,x,y,z,t);
		}
	}
}

///...
void DeepCopy(ScalarField *ScalarField1,ScalarField *ScalarField2,int t)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField2->P(ScalarField1->G(x,y,z),x,y,z,t);
	}	
}


/// Compute the scalar product between the vector fields and put it in ScalarField0 (for which NT=1)
void ScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t,float cste)
{
	int i,x,y,z;
	float temp;
	for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		temp=0.0;
		for (i=0;i<3;i++){temp+=VectorField1->G(i,x,y,z,t)*VectorField2->G(i,x,y,z,t);}
		ScalarField0->P(cste*temp,x,y,z);
	}
}


/// Compute the scalar product between the scalar fields and put it in ScalarField0 (for which NT=1)
void ScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t,float cste)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField0->P(cste * ScalarField1->G(x,y,z,t)*ScalarField2->G(x,y,z,t),x,y,z);
	}	
}


/// Add the scalar product between the vector fields to ScalarField0 (for which NT=1)
void AddScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t)
{
	int i,x,y,z;
	for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		for (i=0;i<3;i++){ScalarField0->Add(VectorField1->G(i,x,y,z,t)*VectorField2->G(i,x,y,z,t),x,y,z);}
	}	
}


/// Add the scalar product between the scalar fields to ScalarField0 (for which NT=1)
void AddScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField0->Add(ScalarField1->G(x,y,z,t)*ScalarField2->G(x,y,z,t),x,y,z);
	}	
}


/// Add  ScalarField1 at time t1 to ScalarField2 at time t2
void AddScalarField(ScalarField *ScalarField1, ScalarField *ScalarField2,float cste, int t1,int t2)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField2->Add(cste*ScalarField1->G(x,y,z,t1),x,y,z,t2);
	}	
}
/// Add  ScalarField1 at time t1 to ScalarField2 at time t2
void AddVectorField(VectorField *VectorField1, VectorField *VectorField2,float cste, int t1,int t2)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		VectorField2->Add(cste*VectorField1->G(i,x,y,z,t1),i,x,y,z,t2);
	}	
}
/// Multiply a vector field by the cste
void MultiplyVectorField(VectorField *VectorField1, float cste,int t)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		VectorField1->P(cste*VectorField1->G(x,y,z,t),i,x,y,z,t);
	}	
}

/// Sum two vector fields and put it in Output.
void SumVectorField(VectorField *VectorField1, VectorField *VectorField2, VectorField *Output, int t1, int t2, int t3, float cste1 ,float cste2)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		Output->P(cste1 * VectorField1->G(i,x,y,z,t1) + cste2 * VectorField2->G(i,x,y,z,t2),i,x,y,z,t3);
	}
}
/// compute the product element by element along each dimension
void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < ScalarField->NZ; z++) for (y = 0; y < ScalarField->NY; y++) for (x = 0; x < ScalarField->NX; x++)
	{
		VectorField2->P(VectorField1->G(i,x,y,z)*ScalarField->G(x,y,z),i,x,y,z);
	}
}


///...
float DotProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, int t1,int t2)
{
	float result=0.0;
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		result += ScalarField1->G(x,y,z,t1) * ScalarField2->G(x,y,z,t2);
	}
	return result;
}
