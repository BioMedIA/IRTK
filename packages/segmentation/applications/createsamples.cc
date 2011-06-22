#include <irtkImage.h>

#ifdef HAS_OPENCV

#include <irtkPointSet.h>

#include <irtkRegistration.h>

#include <cxcore.h>

#include <highgui.h>

#include <nr.h>

/*
reading in input sequence, find mid slice from sequence, using landmark to define interest region
output to jpg files write interest region to positive output info.dat on the end of the file

output 2 images without interest object from the same slice to negtive output and write file name to info.dat
in the corresponding fold.
*/

void usage()
{
	cerr << "Usage: createsamples [inputsequence] [landmark] [positive output directory] [negtive output directory] [objectname]\n" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	
	int i, k, x, y, z, t;
	char buffer[255],posdirectorybuffer[255],negdirectorybuffer[255];
	if(argc < 5){
		usage();
	}

	// Determine how many volumes we have
	irtkGreyImage inputsequence;

	cout << "Reading image: " << argv[1] << endl;
	inputsequence.Read(argv[1]);

	irtkImageAttributes inputSeqAttributes=inputsequence.GetImageAttributes();

	t = inputSeqAttributes._t;

	inputSeqAttributes._t = 1;

	//irtkGreyImage* input = new irtkGreyImage[t];

	// Read first landmark
	irtkPointSet interest;
	//irtkPointSet *rest = new irtkPointSet[t-1];
	cout << "Reading landmark: " << argv[2] << endl;
	interest.ReadVTK(argv[2]);
    //check if landmark's size is 2
	if(interest.Size() != 2){
		cerr<<"Landmark wrong"<<endl;
	    exit(1);
	}
	inputsequence.WorldToImage(interest(0));
	inputsequence.WorldToImage(interest(1));
	int hc_x = (interest(0)._x + interest(1)._x)/2;
	int hc_y = (interest(0)._y + interest(1)._y)/2;

	char* posfilename;
	posfilename = argv[3];

	char* negfilename;
	negfilename = argv[4];

	sprintf(posdirectorybuffer, "%s\\info.dat",posfilename);

	sprintf(negdirectorybuffer, "%s\\info.dat",negfilename);

	char* objectname;
	objectname = argv[5];

	cout << "Write out images" << endl;
	// initialize first image.
    irtkGreyPixel min,max;
	inputsequence.GetMinMax(&min,&max);

	//create jpg image
	IplImage* pimage= cvCreateImage(cvSize(inputSeqAttributes._x, inputSeqAttributes._y), IPL_DEPTH_8U, 1);
	
	for (i = 0; i < 1; i++) {
		cout << "Writing Volume: " << i+1 << " ..." << endl;
		//input[i].Initialize(inputSeqAttributes);
		for ( z = ((inputSeqAttributes._z - 1)/ 2 - round(inputSeqAttributes._z/8.0)); z < ((inputSeqAttributes._z - 1)/ 2 + round(inputSeqAttributes._z/8.0+1)); z++){			

			//write pixel
			for (y = 0; y < inputSeqAttributes._y; y++) {
				for (x = 0; x < inputSeqAttributes._x; x++) {
					//input[i](x, y, z) = inputsequence(x,y,z,i);
					int tmp = (inputsequence(x,y,z,i) * 256 /max);
					pimage->imageData[y*pimage->widthStep + x] = tmp;
				}
			}
			cvEqualizeHist( pimage, pimage );
			// output to positive
			cout << "Output Positive Image " << i << " ..." << endl;
			sprintf(buffer, "%s\\%s%.2d%.2d.jpg",posfilename,objectname,i,z);
			cvSaveImage(buffer,pimage);
			//write to infodata
			ofstream fout(posdirectorybuffer,ios::app);
			sprintf(buffer, "%s%.2d%.2d.jpg",objectname,i,z);
			fout << buffer << " 1 " << round(interest(0)._x) << " " << round(interest(0)._y) 
				<<" " << round(interest(1)._x - interest(0)._x) << " " << round(interest(1)._y - interest(0)._y) << endl;
			fout.close();
			//delete pimage;
			//for (j=0;j<2;j++){
				//output to negtive
				//create jpg image
			IplImage* limg= cvCreateImage(cvSize(hc_x, inputSeqAttributes._y), IPL_DEPTH_8U, 1);
			IplImage* rimg= cvCreateImage(cvSize(inputSeqAttributes._x - hc_x, inputSeqAttributes._y), IPL_DEPTH_8U, 1);
			IplImage* timg= cvCreateImage(cvSize(inputSeqAttributes._x, hc_y), IPL_DEPTH_8U, 1);
			IplImage* bimg= cvCreateImage(cvSize(inputSeqAttributes._x, inputSeqAttributes._y - hc_y), IPL_DEPTH_8U, 1);
			//write pixel
			int ymin[4] = {0,0,0,hc_y};
			int ymax[4] = {inputSeqAttributes._y,inputSeqAttributes._y,hc_y,inputSeqAttributes._y};
			int xmin[4] = {0,hc_x,0,0};
			int xmax[4] = {hc_x,inputSeqAttributes._x,inputSeqAttributes._x,inputSeqAttributes._x};
			for (k=0 ; k<4; k++){
				for (y = ymin[k]; y < ymax[k]; y++) {
					for (x = xmin[k]; x < xmax[k]; x++) {
						int tmp;
						//switch (j){				
						//case 0: tmp = (inputsequence(x,y,z,i) * 256 /max); break;
						//case 1: tmp = (inputsequence(x,y,inputSeqAttributes._z-1,i) * 256 /max); break;
						//}
						tmp = (inputsequence(x,y,z,i) * 256 /max);
						switch(k){
							 case 0: limg->imageData[(y - ymin[k])*limg->widthStep + x-xmin[k]] = tmp; break;
							 case 1: rimg->imageData[(y - ymin[k])*rimg->widthStep + x-xmin[k]] = tmp; break;
							 case 2: timg->imageData[(y - ymin[k])*timg->widthStep + x-xmin[k]] = tmp; break;
							 case 3: bimg->imageData[(y - ymin[k])*bimg->widthStep + x-xmin[k]] = tmp; break;
						}
					}
				}
				switch(k){
					 case 0: cvEqualizeHist( limg, limg ); break;
					 case 1: cvEqualizeHist( rimg, rimg ); break;
					 case 2: cvEqualizeHist( timg, timg ); break;
					 case 3: cvEqualizeHist( bimg, bimg ); break;
				}
				// output to negtive
				cout << "Output Negtive Image " << i <<" " << z << " " << k << " ..." << endl;
				sprintf(buffer, "%s\\%s%.2d%.2d_%d.jpg",negfilename,objectname,i,z,k);
				switch(k){
					case 0: cvSaveImage(buffer,limg); break;
					case 1: cvSaveImage(buffer,rimg); break;
					case 2: cvSaveImage(buffer,timg); break;
					case 3: cvSaveImage(buffer,bimg); break;
				}
				//write to infodata
				ofstream fout(negdirectorybuffer,ios::app);
				//sprintf(buffer, "%s%d%.2d%.2d_%d.jpg",objectname,j,z,i,k);
				sprintf(buffer, "%s%.2d%.2d_%d.jpg",objectname,i,z,k);
				fout << buffer << endl;
				fout.close();
				//delete pimage;
			}
			cvReleaseImage(&limg);
			cvReleaseImage(&rimg);
			cvReleaseImage(&timg);
			cvReleaseImage(&bimg);
		}
	}
	cvReleaseImage(&pimage);
}
#else
void usage()
{
	cerr << "Must be compiled with OpenCV\n" << endl;
	exit(1);
}

int main(int, char **)
{
	usage();
}
#endif
