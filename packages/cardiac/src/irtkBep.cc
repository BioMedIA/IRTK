#include "irtkBep.h"

#ifdef HAS_VTK

void irtkBep::SetOutput(char *outputfilename){
	this->outputfilename=outputfilename;
}

void irtkBep::SetInput(vtkPolyData *bepsurface, vtkPolyData *datasurface){
	this->bepsurface = bepsurface;
	this->datasurface = datasurface;
}

void irtkBep::Initialize(){
	if(bepsurface == NULL ||
		datasurface == NULL || outputfilename == NULL){
			cerr << "please set input and output" << endl;
			exit(1);
	}
}

void irtkBep::Finalize(){
	
}

void irtkBep::EvaluateInfarction (irtkRealImage *threshold, irtkRealImage **input, int number,double timeweight,int cutmode,int connectmode, double regionweight)
{
	double treshold = 0.0001;
	int iterations = 20, n = 2, i;

	irtkEMClassification *classification;
	classification= new irtkEMClassification();
	double *mean, *var, *c;
	double rel_diff;

	//Default settings for Gaussian Mixture parameters
	mean = new double[n];
	var  = new double[n];
	c    = new double[n];

	irtkRealPixel average,std;
	average = input[0]->GetAverage();
	std = input[0]->GetSD();

	mean[0] = average/2.0;
	mean[1] = average + 4*std;
	var[0] = std*std;
	var[1] = 2*std*std;
	for (i=0;i<n;i++)  c[i] = 1.0/(double) n;

	classification->SetInput(*input[0]);
	classification->SetPadding(-1);
	classification->CreateMask();
	classification->InitialiseGMMParameters(n,mean,var,c);

	i=1;
	do {
		cout << "Iteration = " << i << " / " << iterations << endl;
		rel_diff = classification->IterateGMM(i);
		i++;
	} while ((rel_diff>treshold)&&(i<iterations));

	irtkRealImage **tmpatest;
    tmpatest = new irtkRealImage*[2];
    for(i = 0; i < 2; i++){
      tmpatest[i] = new irtkRealImage;
    }
	classification->GetProbMap(0,*tmpatest[0]);
	classification->GetProbMap(1,*tmpatest[1]);

	//add weights for graph cut
	irtkImageGraphCut<float> graphcut;
	graphcut.SetInput(number,input,2,tmpatest);
	graphcut.SetOutput(threshold);
	graphcut.SetMode(cutmode);
	graphcut.Setdt(timeweight);
	graphcut.Run(regionweight,connectmode);
    for(i = 0; i < 2; i++){
        delete tmpatest[i];
    }
    delete []tmpatest;
	delete []mean;
	delete []var;
	delete []c;
	delete classification;
}

void irtkBep::Bullseyeplot(){
  int i,t;
  double bep[17],count[17];

  vtkDoubleArray *beparray = vtkDoubleArray::New();
  vtkDoubleArray *dataarray = vtkDoubleArray::New();
  beparray = (vtkDoubleArray*)bepsurface->GetPointData()->GetScalars();
  dataarray = (vtkDoubleArray*)datasurface->GetPointData()->GetScalars();

  for(i=0; i<17; i++){
	  bep[i] = 0;
	  count[i] = 0;
  }

  vtkPointLocator *pointLocator = vtkPointLocator::New();
  pointLocator->SetDataSet(datasurface);
  pointLocator->BuildLocator();

  for(t=0; t<beparray->GetNumberOfTuples(); t++){
      double point[3];
      vtkIdType ptId;
      bepsurface->GetPoints()->GetPoint (t, point);
      ptId = pointLocator->FindClosestPoint(point);
	  for(i=0; i<17; i++){
		  if(*beparray->GetTuple(ptId) == i){
			  count[i] ++;
			  bep[i] += *dataarray->GetTuple(t);
		  }
	  }
  }

  for(i=0; i<17; i++){
	  if(count[i] == 0)
		  bep[i] = 0;
	  else
		  bep[i] = bep[i]/count[i];
  }
  
  //Output bulleyeplot to file
  ofstream fout(outputfilename,ios::app);
  for(i=0; i<17; i++)
	  fout << bep[i] <<" ";
  fout << endl;
  fout.close();

}

irtkBep::irtkBep(){
	outputfilename = NULL;
	bepsurface = NULL;
	datasurface = NULL;
}

#endif