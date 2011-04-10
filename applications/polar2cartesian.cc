#include <irtkImage.h>
#include <irtkImageToImage.h>
#include <irtkImageFunction.h>
#include <irtkGenericImage.h>

char *input_name = NULL, *output_name = NULL;
#define pi 3.1415926

void usage()
{
	cerr << "Usage: p2c [in] [out] <options>\n";
	cerr << "Where <options> are one or more of the following:\n";
	cerr << "\t<-size x y>      New image dimension \n";
	cerr << "\t<-linear>          Linear interpolation\n";
	cerr << "\t<-bspline>         B-spline interpolation\n";
	cerr << "\t<-cspline>         Cubic spline interpolation\n";
	cerr << "\t<-sinc>            Truncated sinc interpolation\n";
	cerr << "\t<-gaussian sigma>  Gaussian interpolation\n";
	exit(1);
}

int main(int argc, char **argv)
{
	bool ok;
	double xsize, ysize;
	irtkGenericImage<irtkGreyPixel> image,output;
	irtkImageFunction *interpolator = NULL;


	// Check command line
	if (argc < 3) {
		usage();
	}

	// Read input and output names
	input_name  = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;

	// Read input
	image.Read(input_name);

	// Parse remaining parameters
	xsize = 0;
	ysize = 0;
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-size") == 0)) {
			argc--;
			argv++;
			xsize = atof(argv[1]);
			argc--;
			argv++;
			ysize = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-linear") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkLinearInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkBSplineInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-cspline") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkCSplineInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkSincInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-gaussian") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkGaussianInterpolateImageFunction(atof(argv[1]));
			argc--;
			argv++;
			ok = true;
		}
	}
	cout << "Resampling ... "; cout.flush();
	
	

	irtkImageAttributes atr;
	atr = image.GetImageAttributes();
	atr._x=xsize;
	atr._y=ysize;
	output.Initialize(atr);

	//interpolation
	interpolator->SetInput(&image);
	interpolator->Initialize();


	//Resampling Image Center
	double Om=(xsize+1)/2;
	double On=(ysize+1)/2;
	//scale factors
	double sx=(xsize-1)/2;
	double sy=(ysize-1)/2;
	//the ratio range of filled data
	double rMax=1,rMin=0;
	double delR=(rMax-rMin)/(image.GetY()-1);
	double delT=2*pi/image.GetX();

	for(int z=0;z<output.GetZ();z++)
	{
		for(int yi=0;yi<output.GetY();yi++)
		{
			for(int xi=0;xi<output.GetX();xi++)
			{
				double x=(xi-Om)/sx;
				double y=(yi-On)/sy;
				double r=sqrt(x*x+y*y);
				if(r>=rMin&&r<=rMax)
				{
					double t=0;
					if (x>0&&y>=0)
						t=atan(y/x);
					if (x>0&&y<0)
						t=atan(y/x)+2*pi;
					if (x<0)
						t=atan(y/x)+pi;
					if (x==0&&y>=0)
						t=pi/2;
					if (x==0&&y<0)
						t=3*pi/2;;
					if (x==0&&y==0)
						t=0;
					double ri = (r - rMin)/delR;
					double ti =	t/delT;
					double v=0;
					//OCT image
					v=interpolator->Evaluate(ti,ri,z);
					output.PutAsDouble(xi,yi,z,v);
				}
			}
		}
		cout<<z+1<<" ";
	} 



	cout << "done" << endl;

	// Save result
	output.Write(output_name);

	return 0;
}