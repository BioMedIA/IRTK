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

#include <irtkGaussianBlurring.h>

char *cine_name = NULL;
char *out_ED_name = NULL;
char *out_ES_name = NULL;

void usage()
{
	cerr << "Usage: cardiacphasedetection [cine] [ED output] [ES output]" << endl;
	exit(1);
}

int main( int argc, char** argv )
{
	int ok,i,j,k,t,esphase,frames;
    short cine_max,cine_min,cinedis;
	double *similarity,*smoothsimilarity,dif;
	// Check command line
	if (argc < 4) {
		usage();
	}

	// Parse source and target images
	cine_name = argv[1];
	argc--;
	argv++;
	out_ED_name = argv[1];
	argc--;
	argv++;
    out_ES_name = argv[1];
    argc--;
    argv++;

    // Create images
    irtkGreyImage cine(cine_name);
    irtkGreyImage blured;
    blured.Initialize(cine.GetImageAttributes());

    irtkGaussianBlurring<irtkGreyPixel> gaussianBlurring(2);
    gaussianBlurring.SetInput (&cine);
    gaussianBlurring.SetOutput(&blured);
    gaussianBlurring.Run();
    
    irtkImageAttributes atr = cine.GetImageAttributes();
    similarity = new double[atr._t];
    smoothsimilarity = new double[atr._t];
    frames = atr._t;
    atr._t = 1;

    irtkGreyImage out_ED(atr);
    irtkGreyImage out_ES(atr);
    out_ED = cine.GetFrame(0);
    // Create similarity
    blured.GetMinMax(&cine_min,&cine_max);
    cinedis = cine_max - cine_min;
    for(i = 0; i < frames; i++){
        similarity[i] = 0;
        smoothsimilarity[i] = 0;
    }
    // Evaluate similarity
    for ( t = 0; t < frames; t++){
        for ( k = 0; k < cine.GetZ(); k++){
            for ( j = 0; j< cine.GetY(); j++){
                for ( i = 0; i<cine.GetX(); i++){
                   dif = (blured.GetAsDouble(i,j,k,t) - blured.GetAsDouble(i,j,k,0))/cinedis;
                   similarity[t] += dif*dif;
                }
            }
        }
        cout<<"similarity : "<<similarity[t]<<endl;
    }
    for(i = 0; i < frames; i++)
        similarity[i] = sqrt(similarity[i]);
    // Smooth similarity
    for(i = 1; i < frames - 1; i++)
        smoothsimilarity[i] = (similarity[i-1] + similarity[i] + similarity[i+1])/3;
    // Find min similarity
    dif = 0;
    for(i = 0; i < frames; i++){
        if(dif < smoothsimilarity[i]){
            dif = smoothsimilarity[i];
            esphase = i;
        }
    }

    cout<<"ES phase is: "<<esphase<<endl;

    out_ES = cine.GetFrame(esphase);
    // Output
    out_ED.Write(out_ED_name);
    out_ES.Write(out_ES_name);
    delete []smoothsimilarity;
    delete []similarity;
}
