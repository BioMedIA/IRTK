#include <irtkGraphCutSegmentation_4D.h>
#include <irtkImage.h>
#include <string>
#include <irtkHistogram_1D.h>
#include <irtkDilation.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <irtkSubcorticalSegmentation_4D.h>
#include <irtkGradientImage.h>

irtkSubcorticalSegmentation_4D::irtkSubcorticalSegmentation_4D(){
	useMask = false;
}

irtkSubcorticalSegmentation_4D::irtkSubcorticalSegmentation_4D(string name, char * structureFile){
	readStructureNames(structureFile);
	imageName = name;
 	useMask = false;
}

void irtkSubcorticalSegmentation_4D::readParameters(char * filename){
	string line;
  ifstream from(filename);

  if (!from) {
    cerr << "Can't open file " << filename
    << endl;
    exit(1);
  }

  if (from.is_open())
  {
    while (! from.eof() )
    {
      getline (from,line);
      this->Read(line);

    }
    from.close();
  }
}

void irtkSubcorticalSegmentation_4D::readStructureNames(char * filename){
  string line;
  ifstream from(filename);
  if (!from) {
    cerr << "Can't open file " << filename
    << endl;
    exit(1);
  }
  if (from.is_open())
  {
    while (! from.eof() )
    {
    	getline (from,line);
       numTissues++;
    }
    from.close();
  }
  numTissues /= 2;
  names = new string*[numTissues];
  structValues = new int*[numTissues];
  ifstream 	from2(filename);  
  cout << "Read " << numTissues << " Structures:" << endl;
  int counter1 = 0;
  int counter2 = 0;
  int pos = 0;
  if (from2.is_open())
  {
    while (! from2.eof() )
    {
    	structValues[counter2] = new int;
    	names[counter1] = new string;
        getline (from2,line);
       
   		if(pos % 2 == 0){
        	*names[counter1++] = line;
        	pos ++;
        	 cout << *names[counter1-1] << ", ";
   		}
        else{
        	*structValues[counter2++] = atoi(line.c_str());
        	pos++;
        	 cout << *structValues[counter2-1] << " " << 	endl;
        }
    }
    from2.close();
  }

}

void irtkSubcorticalSegmentation_4D::generateGaussians(double minPerc){

    int _padding = -1;
	vector<double> mi_num(numTissues);
	vector<double> denom(numTissues);
	vector<double> sigma_num(numTissues);
	vector<double> mi(numTissues);
	vector<double> *structureValues = new vector<double>[numTissues];
	irtkRealPixel *ptr = input.GetPointerToVoxels();
	int min[numTissues];
	int max[numTissues];
	irtkRealPixel *maxPerc = new irtkRealPixel[numTissues];
	for(int k = 0; k < numTissues; k++){
		min[k]=32767;
	  	max[k]=-32768;
	  	irtkRealPixel tmp;
	  	atlasI[k]->GetMinMax(&tmp,&maxPerc[k]);
	  	cout << "Structure " << k << ": Max atlas value: " << maxPerc[k] << "." << endl;
	}
	irtkRealPixel ** tissuePtr = new irtkRealPixel*[3];
	for(int i = 0; i < 3; i++){
		tissuePtr[i] = new irtkRealPixel;
		tissuePtr[i] = tissuePriors[i]->GetPointerToVoxels();
	}
	for(int i = 0; i < 3 ; i++){
		counter[i]=tissueMy[i]=tissueS[i]=0;
	}
	int atlasVoxelNum = atlasI[0]->GetNumberOfVoxels();
    int inputVoxelNum =  input.GetNumberOfVoxels();

	int whiteVoxelNum = tissuePriors[0]->GetNumberOfVoxels();
	if(whiteVoxelNum != atlasVoxelNum || atlasVoxelNum != inputVoxelNum   || inputVoxelNum != atlasVoxelNum){
		cerr << "Different Voxel Numbers!" << endl;
		exit(1);
 
	} 
	irtkRealPixel **atlasPtr;
	atlasPtr = new irtkRealPixel*[numTissues];
	for(int i = 0; i < numTissues; i++){
		atlasPtr[i] = new irtkRealPixel;
		atlasPtr[i] = atlasI[i]->GetPointerToVoxels();
	}
	for (int i = 0; i < input.GetNumberOfVoxels(); i++){

			for(int ii = 0; ii < 3; ii++){
				counter[ii] += (double)*tissuePtr[ii];
				tissueMy[ii] += *ptr * (double)*tissuePtr[ii];
				tissuePtr[ii]++;
			}
	  		for (int k = 0; k < numTissues; k++){
				if(*atlasPtr[k] >= maxPerc[k]*minPerc && *ptr > _padding){
					structureValues[k].push_back(*ptr);
					if(*ptr>max[k])
						max[k] = *ptr;
					if(*ptr<min[k])
						min[k] = *ptr;
				}
				atlasPtr[k]++;
	  		}
	    ptr++;
	}

	for(int i = 0; i < 3; i++){
		tissueMy[i] /= counter[i];
	}

	ptr = input.GetPointerToVoxels();
	for(int ii = 0; ii < 3; ii++){
		tissuePtr[ii] = tissuePriors[ii]->GetPointerToVoxels();
	}
	for (int i = 0; i < input.GetNumberOfVoxels(); i++){
		for(int ii = 0; ii < 3; ii++){
			if(*tissuePtr[ii] > 0){
				tissueS[ii] += (tissueMy[ii]-*ptr)*(tissueMy[ii]-*ptr);
			}
			tissuePtr[ii]++;
		}
	ptr++;
	}
	cout << endl << "Tissue Classes: " << endl; 
	for(int i = 0; i < 3 ; i++){
		cout   << i+1 << ": my: " << tissueMy[i] << " sigma: " << tissueS[i] /counter[i]<< endl;
	}



	for (int k = 0; k < numTissues; k++){
		for(unsigned int i = 0; i < structureValues[k].size(); i++){
			double val = structureValues[k][i];
				mi_num[k] += val;
				denom[k]  += 1;
		}
	    if (denom[k] != 0){
	      	mi[k] = mi_num[k] / denom[k];
    }
    else{
      	cerr << "Division by zero while computing tissue mean!" << endl;
      	exit(1);
    	}
  	}
	    for (int k = 0; k <numTissues; k++){
	    	for(unsigned int i = 0; i < structureValues[k].size(); i++){
	    		double val = structureValues[k][i];
	    			sigma_num[k] += ((val - mi[k]) * (val - mi[k]));
	    	}

	   }
 
    cout << endl << "Structures: " << endl ;
	for (int k = 0; k <numTissues; k++){
		_my[k] = mi[k];
		_sigma[k] = sigma_num[k] / denom[k];
		cout << *names[k] << ": my: " << k << " " << _my[k] << " sigma: " << _sigma[k] << "...no of voxels considered: " << denom[k]<< endl;
		cout << *names[k] << ": my: " << k << " " << _my[k] << " sigma: " << _sigma[k] << "...no of voxels considered: " << denom[k]<< endl;
  	}
  	cout << endl;

  	delete []structureValues;
  	delete []maxPerc;
}

void irtkSubcorticalSegmentation_4D::doSegmentation(double l, double g, double c, double sigmaG){
	cout << "Do segmentation with c = " << c << ", lambda = " << l << ", alpha = " << g << ", sigmaG = " << sigmaG << endl;
	cout << "Get Gradient Magnitude " << endl;
	irtkGradientImageX<irtkRealPixel> gradientX;
	irtkGradientImageY<irtkRealPixel> gradientY;
	irtkGradientImageZ<irtkRealPixel> gradientZ;
	int _xDim = input.GetX();
	int _yDim = input.GetY();
	int _zDim = input.GetZ();
	irtkGenericImage<irtkRealPixel> *_xGrad;
	irtkGenericImage<irtkRealPixel> *_yGrad;
	irtkGenericImage<irtkRealPixel> *_zGrad;
	int _tDim = input.GetT();
 	_xGrad = new irtkGenericImage<irtkRealPixel>(_xDim, _yDim, _zDim, _tDim);
	_yGrad = new irtkGenericImage<irtkRealPixel>(_xDim, _yDim, _zDim, _tDim);
	_zGrad = new irtkGenericImage<irtkRealPixel>(_xDim, _yDim, _zDim, _tDim);
	gradientX.SetInput(&input);
	gradientY.SetInput(&input);
	gradientZ.SetInput(&input);
	gradientX.SetOutput(_xGrad);
	gradientY.SetOutput(_yGrad);
	gradientZ.SetOutput(_zGrad);
	gradientX.Run();
	gradientY.Run();
	gradientZ.Run();
	irtkGenericImage<irtkRealPixel> gradMagnitude;
	gradMagnitude.Initialize(input.GetImageAttributes());
        
  	for(int k = 0; k < 3; k++){
		cout << k << " " << tissueMy[k] << " " << tissueS[k] << " " << counter[k] << endl;
	}	
	irtkRealImage **seg_results;
	seg_results = new irtkRealImage*[numTissues];

	MOG mog(3, 0);

	for(int i = 0; i < 3; i++){
		mog.setTissueProb(i, tissuePriors[i], tissueMy[i], tissueS[i]/counter[i]);
	}
	for(int j = 0; j < numTissues; j++){
		irtkGraphCutSegmentation_4D *gcs;
		cout << endl;
		cout << j << " Preprocessings " << endl;
	 	irtkRealImage temp;
	 	temp.Initialize(input.GetImageAttributes());
 		irtkRealPixel * tempPtr = temp.GetPointerToVoxels();
 		irtkRealPixel * xPtr = _xGrad->GetPointerToVoxels();
		irtkRealPixel * yPtr = _yGrad->GetPointerToVoxels();
		irtkRealPixel * zPtr = _zGrad->GetPointerToVoxels();
		irtkRealPixel * gradPtr = gradMagnitude.GetPointerToVoxels();
		irtkRealPixel * atlasPtr = atlasI[j]->GetPointerToVoxels();
		for(int i = 0; i < input.GetNumberOfVoxels(); i++){ 
				double gradMagnVal = sqrt(*xPtr**xPtr+*yPtr**yPtr+*zPtr**zPtr);
				*gradPtr = gradMagnVal;
				if(*atlasPtr > 0){
					*tempPtr = 1;
				}
				tempPtr++;
				xPtr++;
				yPtr++;
				zPtr++;
				gradPtr++;
				atlasPtr++;	
		}
	 	irtkDilation<irtkRealPixel> dilation;
	  	dilation.SetInput(&temp);
	  	dilation.SetOutput(&temp);
	 	dilation.Run();
	 	irtkRealPixel *ptr = input.GetPointerToVoxels();
	 	irtkRealPixel *ptr2 = temp.GetPointerToVoxels();
	 	for(int i = 0; i < input.GetNumberOfVoxels(); i++){
	 		if(*ptr2 == 1)
	 			*ptr2 = *ptr;
	 		else
	 			*ptr2 = -1;
	 		ptr++;
	 		ptr2++;
	 	}
		gcs = new irtkGraphCutSegmentation_4D(2, *atlasI[j]);
		if(useMask){
			cout << "Using mask" << endl;
			gcs->setMask(mask);
		}
		gcs->GenerateGaussians(0, _my[j], _sigma[j]);
		cout << "j: " << j << "My: " << _my[j] << ". Sigma: " << _sigma[j]  << endl;
		gcs->SetMog(mog);
		gcs->setParameters(l, g, c, sigmaG);
		gcs->SetTissuePriors(tissuePriors);
		gcs->SetInput(temp, gradMagnitude);
		gcs->Iterate();
		string segName = outputDir + "/" + *names[j] + ".nii.gz";
		seg_results[j] = new irtkRealImage;
		*seg_results[j] = gcs->GetSegmentation();
		delete gcs;
	}

	cout << endl << "Fusing individual segmentations" << endl;
	irtkRealImage segmentation;
	irtkGreyImage ** output;
	output = new irtkGreyImage*[input.GetT()];
	irtkImageAttributes attributes = input.GetImageAttributes();
	attributes._t = 1;
	for(int i = 0; i < input.GetT(); i++){
		output[i] = new irtkGreyImage();
		
		output[i]->Initialize(attributes);
	}
	segmentation.Initialize(attributes);

	for(int x = 0; x < input.GetX(); x++){
		for(int y = 0; y < input.GetY(); y++){
			for(int z = 0; z < input.GetZ(); z++){
				for(int t = 0; t < input.GetT(); t++){
					for(int i = 0; i < numTissues; i++){
						if(seg_results[i]->Get(x, y, z, t) == 1){						
							output[t]->Put(x, y, z, 0, *structValues[i]);
						}		
					}
				}
			}

		}
	}

	cout << numTissues << endl;
	for(int i = 0; i < numTissues; i++){
		cout << *structValues[i] << endl;

	}
 	cout << "Write output" << endl;
 	for(int i = 0; i < input.GetT(); i++){
 		ostringstream stm;
		stm << i;
		
		string segName = outputDir + "/" + imageName + "_" + stm.str() + ".nii.gz";
		if(input.GetT() == 1){
			segName = outputDir + "/" + imageName + ".nii.gz";
		}
		output[i]->Write(segName.c_str());
 	}
}

bool irtkSubcorticalSegmentation_4D::Read(string line)
{
  bool ok = false;

  if (line.find( "Image Directory") != string::npos) {
  	 int pos = line.find("=");
    string part = line.substr(pos+2, line.length()-pos-3);
 	imageDir = part;
    cout << "Image Directory: " << imageDir<< endl;
    ok = true;
  }

 if (line.find( "Atlas Directory") != string::npos) {
    int pos = line.find("=");
    atlasDir = line.substr(pos+2, line.length()-pos-3);
    cout << "Atlas Directory: " << atlasDir << endl;
    ok = true;
  }

  if (line.find( "Output Directory") != string::npos) {
    int pos = line.find("=");
   outputDir = line.substr(pos+2, line.length()-pos-3);
    cout << "Output Directory: " << outputDir << endl;
    ok = true;
  }



  if (line.find( "Tissue Apriori Directory") != string::npos) {
    int pos = line.find("=");
    tissuePriorDir= line.substr(pos+2, line.length()-pos-3);
    cout << "Tissue Apriori Directory: " << tissuePriorDir << endl;
    ok = true;
  }

  return ok;
}

void irtkSubcorticalSegmentation_4D::setMask(const char * fileName){
	mask.Read(fileName);
	cout << "Read mask " << fileName << endl;
	useMask = true;
}

void irtkSubcorticalSegmentation_4D::init(irtkRealImage input_image, irtkRealImage ** structure_atlases, irtkRealImage ** tissue_priors, int num_tissues, string out_dir, string ** input_names, int ** struct_values, string image_name){
	cout << "Init via command line arguments" << endl;
	input = input_image;
	atlasI = structure_atlases;
	tissuePriors = tissue_priors;
	numTissues = num_tissues;
	outputDir = out_dir;
	names = input_names;
	structValues = struct_values;
	imageName = image_name;
	irtkRealPixel t_min, t_max;
	tissuePriors[0]->GetMinMax(&t_min, &t_max);
	if(t_max>1){
		for(int i = 0; i < 3; i++){
			irtkRealPixel *ptr = tissuePriors[i]->GetPointerToVoxels();
			for(int j = 0; j < tissuePriors[i]->GetNumberOfVoxels(); j++){
				*ptr /= 255;
				ptr++;
			}
		}
	}

}

void irtkSubcorticalSegmentation_4D::init(){
	cout << "Init via parameter file" << endl;
	string inputName = imageDir + imageName + ".nii.gz";
	input.Read(inputName.c_str());
	atlasI = new irtkRealImage*[numTissues];
	for(int i = 0; i < numTissues; i++){
		atlasI[i] = new irtkRealImage;
		string atlasName = atlasDir + "/" +imageName + "_" + *names[i] + ".nii.gz";
		atlasI[i]->Read(atlasName.c_str());		
	}
	string wmName = tissuePriorDir + "/" + imageName+ "_white.nii.gz";
	string csfName = tissuePriorDir + "/"  + imageName+ "_csf.nii.gz";
	string greyName = tissuePriorDir + "/"  + imageName+ "_grey.nii.gz";
	tissuePriors = new irtkRealImage*[3];
	tissuePriors[0] = new irtkRealImage;
	tissuePriors[1] = new irtkRealImage;
	tissuePriors[2] = new irtkRealImage;

	tissuePriors[0]->Read(wmName.c_str());
	tissuePriors[1]->Read(csfName.c_str());
	tissuePriors[2]->Read(greyName.c_str());
	irtkRealPixel t_min, t_max;
        tissuePriors[0]->GetMinMax(&t_min, &t_max);
        if(t_max>1){	
		for(int i = 0 ; i < 3; i++){
			irtkRealPixel *ptr = tissuePriors[i]->GetPointerToVoxels();
			for(int j = 0; j < tissuePriors[i]->GetNumberOfVoxels(); j++){
				*ptr /= 255;
				ptr++;
			}
		}
	}

}
