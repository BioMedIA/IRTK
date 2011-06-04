#include "irtkTagFunction.h"

double tox, toy, toz, sox, soy, soz;
irtkMatrix tmat(4, 4);
irtkMatrix smat(4, 4);
irtkMatrix tempMat, transfMat;

irtkTagFunction::irtkTagFunction(void)
{
	this->_interest = NULL;
	this->_source = NULL;
	this->_target = NULL;
	this->_threshold = NULL;
	this->_transformation = NULL;
	this->_registration = NULL;
	this->_interestaffine = NULL;
	this->_maxmap = NULL;
	this->sampled = 0;
	this->_exemption = NULL;
}

irtkTagFunction::~irtkTagFunction(void)
{
}
// SetInput of the function
inline void irtkTagFunction::SetInput(irtkGreyImage* interest, irtkRealImage* threshold, irtkGreyImage* target, irtkGreyImage* source, irtkAffineTransformation* interestaffine)
{
	this->_interest = interest;
	this->_source = source;
	this->_target = target;
	this->_threshold = threshold;
	this->_interestaffine = interestaffine;
}

// SetOutput pointset
inline irtkPointSet irtkTagFunction::GetOutput()
{
	return this->_p2;
}

inline irtkPointSet irtkTagFunction::GetPointSet(void){
	return this->_p1;
}

void irtkTagFunction::InitializeTransformation(irtkGreyImage &target, irtkGreyImage &source)
{
	tox = toy = toz = sox = soy = soz = 0.0;
	tmat.Ident();
	smat.Ident();

    this->_transformation = new irtkAffineTransformation;

	// Registration
	this->_registration = new irtkImageAffineRegistration;

	// Set input and output for the registration filter
	this->_registration->SetOutput(this->_transformation);

	cout << "Centering ... ";
    // Place the voxel centre at the world origin.
    target.GetOrigin(tox, toy, toz);
    source.GetOrigin(sox, soy, soz);
    target.PutOrigin(0.0, 0.0, 0.0);
    source.PutOrigin(0.0, 0.0, 0.0);

    // Update the transformation accordingly.
    tmat(0, 3) = tox;
    tmat(1, 3) = toy;
    tmat(2, 3) = toz;
    smat(0, 3) = -1.0 * sox;
    smat(1, 3) = -1.0 * soy;
    smat(2, 3) = -1.0 * soz;

    transfMat = _transformation->GetMatrix();
    tempMat   = transfMat * tmat;
    tempMat   = smat * tempMat;

    _transformation->PutMatrix(tempMat);
    cout << "done" << endl;

}

void irtkTagFunction::FinalizeTransformation(irtkGreyImage &target, irtkGreyImage &source)
{

	target.PutOrigin(tox, toy, toz);
	source.PutOrigin(sox, soy, soz);

	tmat(0, 3) = -1.0 * tox;
	tmat(1, 3) = -1.0 * toy;
	tmat(2, 3) = -1.0 * toz;
	smat(0, 3) = sox;
	smat(1, 3) = soy;
	smat(2, 3) = soz;

	transfMat = _transformation->GetMatrix();
	tempMat   = transfMat * tmat;
	tempMat   = smat * tempMat;

	_transformation->PutMatrix(tempMat);

}


void irtkTagFunction::Initialize(void)
{
	if(this->sampled == 0){
		this->Sample();
	    this->_p2 = this->_p1;
	}else
		this->_p2 = this->_p1;
}

void irtkTagFunction::Finalize(void)
{
	delete this->_transformation;
	delete this->_registration;
}

void irtkTagFunction::Sample(void)
{
	if(this->_interest == NULL || this->_threshold == NULL || this->_interestaffine == NULL){
		cerr<<"interest region not given can't sample tags"<<endl;
	    exit(1);
	}
    //define pointers
	irtkGreyPixel *ptr2target,*ptr2maxmap, *ptr2exemption;
	irtkRealPixel *ptr2threshold;
	int i,j,k,t;
	int i1,j1;
	int times, average;
	//caculate interest region;
	irtkPoint p1,p2;
	p1._x = this->_interestaffine->Get(0) - this->_interestaffine->Get(6)/this->_target->GetXSize()/2;
    p1._y = this->_interestaffine->Get(1) - this->_interestaffine->Get(7)/this->_target->GetYSize()/2;
    p1._z = 0;

    p2._x = this->_interestaffine->Get(0) + this->_interestaffine->Get(6)/this->_target->GetXSize()/2;
    p2._y = this->_interestaffine->Get(1) + this->_interestaffine->Get(7)/this->_target->GetYSize()/2;
    p2._z = _target->GetZ();

    //before sample find out the tag spacing using fourier transformation	
	/*
	irtkGreyPixel min,max;
	int x,y,z;
	//create opencv images from _target
	irtkGreyImage fimg(this->_target->GetRegion(1,1,0,this->_target->GetX()-1,this->_target->GetY()-1,this->_target->GetZ()));
	//= this->_target->GetRegion(p1._x,p1._y,p1._z,p2._x,p2._y,p2._z);
	fimg.GetMinMax(&min,&max);

	IplImage* im= cvCreateImage(cvSize(fimg.GetX(), fimg.GetY()), IPL_DEPTH_8U, 1);

	//write pixel
	z = (fimg.GetZ() - 1)/ 2;
	for (y = 0; y < fimg.GetY(); y++) {
		for (x = 0; x < fimg.GetX(); x++) {
			int tmp = (fimg(x,y,z,0) * 256 /max);
			im->imageData[y*im->widthStep + x] = tmp;
		}
	}

	//cvEqualizeHist( im, im );

	//do fourier transformation
	IplImage * realInput;
    IplImage * imaginaryInput;
    IplImage * complexInput;
    int dft_M, dft_N;
    CvMat* dft_A, tmp;
    IplImage * image_Re;
    IplImage * image_Im;
    double m, M;

	realInput = cvCreateImage( cvGetSize(im), IPL_DEPTH_64F, 1);
    imaginaryInput = cvCreateImage( cvGetSize(im), IPL_DEPTH_64F, 1);
    complexInput = cvCreateImage( cvGetSize(im), IPL_DEPTH_64F, 2);
 
    cvScale(im, realInput, 1.0, 0.0);
    cvZero(imaginaryInput);
    cvMerge(realInput, imaginaryInput, NULL, NULL, complexInput);
 
    dft_M = cvGetOptimalDFTSize( im->height - 1 );
    dft_N = cvGetOptimalDFTSize( im->width - 1 );
 
    dft_A = cvCreateMat( dft_M, dft_N, CV_64FC2 );
    image_Re = cvCreateImage( cvSize(dft_N, dft_M), IPL_DEPTH_64F, 1);
    image_Im = cvCreateImage( cvSize(dft_N, dft_M), IPL_DEPTH_64F, 1);
 
    // copy A to dft_A and pad dft_A with zeros
    cvGetSubRect( dft_A, &tmp, cvRect(0,0, im->width, im->height));
    cvCopy( complexInput, &tmp, NULL );
    if( dft_A->cols > im->width )
    {
        cvGetSubRect( dft_A, &tmp, cvRect(im->width,0, dft_A->cols - im->width, im->height));
        cvZero( &tmp );
    }
 
    // no need to pad bottom part of dft_A with zeros because of
    // use nonzero_rows parameter in cvDFT() call below
 
    cvDFT( dft_A, dft_A, CV_DXT_FORWARD, complexInput->height );
 
    cvNamedWindow("win", 0);
    cvNamedWindow("magnitude", 0);
	cvNamedWindow("imaginary", 0);
    cvShowImage("win", im);
 
    // Split Fourier in real and imaginary parts
    cvSplit( dft_A, image_Re, image_Im, 0, 0 );
 
    // Compute the magnitude of the spectrum Mag = sqrt(Re^2 + Im^2)
    cvPow( image_Re, image_Re, 2.0);
    cvPow( image_Im, image_Im, 2.0);
    cvAdd( image_Re, image_Im, image_Re, NULL);
    cvPow( image_Re, image_Re, 0.5 );
 
    // Compute log(1 + Mag)
    cvAddS( image_Re, cvScalarAll(1.0), image_Re, NULL ); // 1 + Mag
    cvLog( image_Re, image_Re ); // log(1 + Mag)
 
 
    // Rearrange the quadrants of Fourier image so that the origin is at
    // the image center
	cvShiftDFT( image_Im, image_Im);
    cvShiftDFT( image_Re, image_Re );
 
    cvMinMaxLoc(image_Re, &m, &M, NULL, NULL, NULL);
    cvScale(image_Re, image_Re, 1.0/(M-m), 1.0*(-m)/(M-m));
    cvShowImage("magnitude", image_Re);
	cvShowImage("imaginary", image_Im);
 
    cvWaitKey(-1);
	//find out the frequency of the tag spacing
	//the above code is test only now */

	//prepare max map image
	this->_maxmap = new irtkGreyImage(*this->_target);
	ptr2threshold = this->_threshold->GetPointerToVoxels();
	ptr2target = this->_target->GetPointerToVoxels();
	ptr2maxmap = this->_maxmap->GetPointerToVoxels();
	for (t=0; t<this->_target->GetT();t++){
		for (k = 0; k < this->_target->GetZ(); k++) {
			for (j = 0; j < this->_target->GetY(); j++) {
				for (i = 0; i <this->_target->GetX(); i++) {
					times = *ptr2threshold - 2;
					if(times > 0)
						*ptr2maxmap = times*(*ptr2target);
					else
						*ptr2maxmap = 0;
					ptr2threshold++;
					ptr2target++;
					ptr2maxmap++;
				}
			}
		}
	}
	//debug
	//this->_maxmap->Write("maxmap.gipl");
	//create exception image
	this->_exemption = new irtkGreyImage(this->_target->GetImageAttributes());
	average = this->_maxmap->GetAverage();

	//loop over image
	for(t=0; t<this->_target->GetT(); t++){
		for (k = round(p1._z)+1; k < round(p2._z)-1; k++) {
			for (j = round(p1._y); j < round(p2._y); j++) {
				for (i = round(p1._x); i < round(p2._x); i++) {
					ptr2threshold = this->_threshold->GetPointerToVoxels(i,j,k,t);
					ptr2exemption = this->_exemption->GetPointerToVoxels(i,j,k,t);
					
					// check if meet the tag profile
					if(this->_maxmap->GetAsDouble(i,j,k,t) > average 
						&& this->_maxmap->GetAsDouble(i-1,j,k,t) > average
						&& this->_maxmap->GetAsDouble(i+1,j,k,t) > average
						&& this->_maxmap->GetAsDouble(i,j-1,k,t) > average
						&& this->_maxmap->GetAsDouble(i,j+1,k,t) > average
						&& *ptr2threshold > 2 && *ptr2exemption == 0){
						irtkPoint tmp(i,j,k);
						this->_maxmap->ImageToWorld(tmp);
						this->_maxmap->GetMaxPosition(tmp,2,t);
						this->_maxmap->GravityCenter(tmp,2,t);

						//exempt near by pixels
						this->_p1.Add(tmp);
						irtkPoint extmp(tmp);
						this->_maxmap->WorldToImage(extmp);
						for (j1 = round(extmp._y) - 2; j1 <  round(extmp._y) + 3; j1++){
							for (i1 = round(extmp._x) - 2; i1 <  round(extmp._x) + 3; i1++){
								this->_exemption->PutAsDouble(i1,j1,k,t,10);
							}
						}
					}
				}
			}
		}
	}

	if(this->_p1.Size() != 0 )
	  this->sampled = 1;

	delete this->_maxmap;
	delete this->_exemption;
}

void irtkTagFunction::Track(int toggle, int ds)
{
	int j;
	this->Initialize();

	if (sampled == 0){
		cerr<<"tag detection failed"<<endl;
        exit(1);
	}

	for (j = 0; j < this->_p2.Size(); j++){
		this->_target->WorldToImage(this->_p2(j));

		int tx1,tx2,ty1,ty2,tz1,tz2;
		tx1 = round(this->_p2(j)._x - 5); tx2 = round(this->_p2(j)._x + 6);
		ty1 = round(this->_p2(j)._y - 5); ty2 = round(this->_p2(j)._y + 6);
		tz1 = round(this->_p2(j)._z - 5); tz2 = round(this->_p2(j)._z + 6);
		if(tx1 < 0) tx1 = 0; if(tx1 > _target->GetX()-1) tx1 = _target->GetX()-1;
		if(tx2 < 1) tx2 = 1; if(tx2 > _target->GetX()) tx2 = _target->GetX();
		if(ty1 < 0) ty1 = 0; if(ty1 > _target->GetY()-1) ty1 = _target->GetY()-1;
		if(ty2 < 1) ty2 = 1; if(ty2 > _target->GetY()) ty2 = _target->GetY();
		if(tz1 < 0) tz1 = 0; if(tz1 > _target->GetZ() - 1) tz2 = _target->GetZ() - 1;
		if(tz2 < 1) tz2 = 1; if(tz2 > _target->GetZ()) tz2 = _target->GetZ();
		irtkGreyImage target(this->_target->GetRegion(tx1,ty1,tz1,tx2,ty2,tz2));

		//if two image is the same we don't need this
		if (this->_target->GetWorldToImageMatrix() != this->_source->GetWorldToImageMatrix()){
		  this->_target->ImageToWorld(this->_p2(j));
		  this->_source->WorldToImage(this->_p2(j));
		}

		tx1 = round(this->_p2(j)._x - 8); tx2 = round(this->_p2(j)._x + 9);
		ty1 = round(this->_p2(j)._y - 8); ty2 = round(this->_p2(j)._y + 9);
		tz1 = round(this->_p2(j)._z - 8); tz2 = round(this->_p2(j)._z + 9);
		if(tx1 < 0) tx1 = 0; if(tx1 > _target->GetX()-1) tx1 = _target->GetX()-1;
		if(tx2 < 1) tx2 = 1; if(tx2 > _target->GetX()) tx2 = _target->GetX();
		if(ty1 < 0) ty1 = 0; if(ty1 > _target->GetY()-1) ty1 = _target->GetY()-1;
		if(ty2 < 1) ty2 = 1; if(ty2 > _target->GetY()) ty2 = _target->GetY();
		if(tz1 < 0) tz1 = 0; if(tz1 > _target->GetZ() - 1) tz2 = _target->GetZ() - 1;
		if(tz2 < 1) tz2 = 1; if(tz2 > _target->GetZ()) tz2 = _target->GetZ();
		
		irtkGreyImage source(this->_source->GetRegion(tx1,ty1,tz1,tx2,ty2,tz2));
		//target.Write("target.gipl");
		//source.Write("source.gipl");
		//cleaar target source image attributes
		irtkImageAttributes pt,ps;
		pt._x = target.GetX(); pt._y = target.GetY(); pt._z = target.GetZ();
		ps._x = source.GetX(); ps._y = source.GetY(); ps._z = source.GetZ();
		pt._xorigin = (source.GetX()-1.0)/2.0; pt._yorigin = (source.GetY()-1.0)/2.0; pt._zorigin = (source.GetZ()-1.0)/2.0;
        ps._xorigin = (source.GetX()-1.0)/2.0; ps._yorigin = (source.GetY()-1.0)/2.0; ps._zorigin = (source.GetZ()-1.0)/2.0;
		  // Calculate origin
		target.PutOrientation(pt._xaxis,pt._yaxis,pt._zaxis);
		target.PutOrigin(pt._xorigin,pt._yorigin,pt._zorigin);
		target.PutPixelSize(pt._dx,pt._dy,pt._dz,pt._dt);
		source.PutOrientation(ps._xaxis,ps._yaxis,ps._zaxis);
		source.PutOrigin(ps._xorigin,ps._yorigin,ps._zorigin);
		source.PutPixelSize(ps._dx,ps._dy,ps._dz,ps._dt);
		//target.Write("target.gipl");
		//source.Write("source.gipl");
		//Initialize Transformation
		this->InitializeTransformation(target,source);
		//set input
		this->_registration->SetInput(&target, &source);
		this->_registration->GuessParameter();
		//this->_transformation->Reset();
		//run this->_registration without output.
		this->_registration->Run();
		this->FinalizeTransformation(target,source);
		//transforme in image space
		//write interestregion
		//irtkCofstream to;
		//to.Open("affinetransformation.dof");
		//this->_transformation->Write(to);
		irtkPoint tmp(pt._xorigin,pt._yorigin,pt._zorigin)
			,tmp2(pt._xorigin,pt._yorigin,pt._zorigin);
		this->_transformation->Transform(tmp._x,tmp._y,tmp._z);
		this->_p2(j) = tmp - tmp2 + this->_p2(j);
		this->_source->ImageToWorld(this->_p2(j));

		// center the tag in window if asked.
		if(toggle == 2){
			this->_source->GravityCenter(this->_p2(j),ds);
		}

	}

	this->Finalize();
}

inline void irtkTagFunction::SetPointSet(irtkPointSet& p1)
{
	this->_p1 = p1;
	this->sampled = 1;
}
#ifdef HAS_OPENCV
void cvShiftDFT(CvArr * src_arr, CvArr * dst_arr )
{
    CvMat * tmp;
    CvMat q1stub, q2stub;
    CvMat q3stub, q4stub;
    CvMat d1stub, d2stub;
    CvMat d3stub, d4stub;
    CvMat * q1, * q2, * q3, * q4;
    CvMat * d1, * d2, * d3, * d4;
 
    CvSize size = cvGetSize(src_arr);
    CvSize dst_size = cvGetSize(dst_arr);
    int cx, cy;
 
    if(dst_size.width != size.width || 
       dst_size.height != size.height){
        cvError( CV_StsUnmatchedSizes, "cvShiftDFT", "Source and Destination arrays must have equal sizes", __FILE__, __LINE__ );   
    }
 
    if(src_arr==dst_arr){
        tmp = cvCreateMat(size.height/2, size.width/2, cvGetElemType(src_arr));
    }
 
    cx = size.width/2;
    cy = size.height/2; // image center
 
    q1 = cvGetSubRect( src_arr, &q1stub, cvRect(0,0,cx, cy) );
    q2 = cvGetSubRect( src_arr, &q2stub, cvRect(cx,0,cx,cy) );
    q3 = cvGetSubRect( src_arr, &q3stub, cvRect(cx,cy,cx,cy) );
    q4 = cvGetSubRect( src_arr, &q4stub, cvRect(0,cy,cx,cy) );
    d1 = cvGetSubRect( dst_arr, &d1stub, cvRect(0,0,cx,cy) );
    d2 = cvGetSubRect( dst_arr, &d2stub, cvRect(cx,0,cx,cy) );
    d3 = cvGetSubRect( dst_arr, &d3stub, cvRect(cx,cy,cx,cy) );
    d4 = cvGetSubRect( dst_arr, &d4stub, cvRect(0,cy,cx,cy) );
 
    if(src_arr!=dst_arr){
        if( !CV_ARE_TYPES_EQ( q1, d1 )){
            cvError( CV_StsUnmatchedFormats, "cvShiftDFT", "Source and Destination arrays must have the same format", __FILE__, __LINE__ ); 
        }
        cvCopy(q3, d1, 0);
        cvCopy(q4, d2, 0);
        cvCopy(q1, d3, 0);
        cvCopy(q2, d4, 0);
    }
    else{
        cvCopy(q3, tmp, 0);
        cvCopy(q1, q3, 0);
        cvCopy(tmp, q1, 0);
        cvCopy(q4, tmp, 0);
        cvCopy(q2, q4, 0);
        cvCopy(tmp, q2, 0);
    }
}
#endif
