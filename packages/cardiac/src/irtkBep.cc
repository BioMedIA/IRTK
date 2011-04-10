#include "irtkBep.h"

void irtkBep::SetInput(irtkGreyImage& target, irtkGreyImage& source, irtkRealImage& threshold, char* osequence, char* lsequence){
	this->target = target; this->source = source; this->threshold = threshold;
	this->osequence = osequence; this->lsequence = lsequence;
}

void irtkBep::SetInput(irtkGreyImage& target, irtkRealImage& threshold, char* osequence){
	this->target = target; 
	this->threshold = threshold;
	this->osequence = osequence;
}

void irtkBep::SetLandmarks(irtkPointSet& olandmarks, int tip, int aa, int am, int mb, int bottom){
	this->olandmarks = olandmarks; this->tip = tip;
	this->aa=aa; this->am=am; 
	this->mb=mb; this->bottom=bottom;
}

void irtkBep::SetOutput(char *outputfilename,char *maxfilename){
	this->outputfilename=outputfilename;
	this->maxfilename = maxfilename;
}

void irtkBep::SetCompare(irtkRealImage& compare){
	this->compare = compare;
}

void irtkBep::SetSegmentation(irtkRealImage& segmentation){
	this->segmentation = segmentation;
}

void irtkBep::Initialize(){
	double R1,R2,R;
	double K,Cx,Cy;
	double kp,x1,y1,x2,y2;
	int i;
	double* p = new double[3];
	p[0] = 0; p[1] = 0; p[2] = 5;
	//regenerate 0
	for(i=0;i<11;i++)
	  landmarks.Add(p);
	for(i=0;i<16;i++)
	  centerpoints.Add(p);
	delete p;
	landmarks(0)._x = olandmarks(0)._x; landmarks(0)._y = olandmarks(0)._y; 
	//evaluate R
	R1 = sqrt(pow(olandmarks(1)._x - olandmarks(0)._x,2) + pow(olandmarks(1)._y - olandmarks(0)._y,2));
	R2 = sqrt(pow(olandmarks(2)._x - olandmarks(0)._x,2) + pow(olandmarks(2)._y - olandmarks(0)._y,2));
	R = round(R1>R2?R1:R2);
	Cx = olandmarks(0)._x; Cy = olandmarks(0)._y;
	//regenerate 6 9
	K = (olandmarks(1)._y - olandmarks(0)._y)/(olandmarks(1)._x - olandmarks(0)._x);
    solveequation( x1,  x2,  y1,  y2,  Cx,  R,  K,  Cy);

	kp = (Cx-olandmarks(2)._x)*(y1-olandmarks(2)._y) - (Cy-olandmarks(2)._y)*(x1-olandmarks(2)._x);
	if(kp > 0){
		landmarks(6)._x = round(x2); landmarks(6)._y = round(y2);
        landmarks(9)._x = round(x1); landmarks(9)._y = round(y1);
		//landmarks(2)._x = round(x2); landmarks(2)._y = round(y2);
        //landmarks(4)._x = round(x1); landmarks(4)._y = round(y1);
	}else{
		landmarks(6)._x = round(x1); landmarks(6)._y = round(y1);
		landmarks(9)._x = round(x2); landmarks(9)._y = round(y2);
		//landmarks(2)._x = round(x1); landmarks(2)._y = round(y1);
		//landmarks(4)._x = round(x2); landmarks(4)._y = round(y2);
	}

	//regenerate 5 8
	K = (olandmarks(2)._y - olandmarks(0)._y)/(olandmarks(2)._x - olandmarks(0)._x);
    solveequation( x1,  x2,  y1,  y2,  Cx,  R,  K,  Cy);

	kp = (Cx-olandmarks(1)._x)*(y1-olandmarks(1)._y) - (Cy-olandmarks(1)._y)*(x1-olandmarks(1)._x);
	if(kp > 0){
		landmarks(5)._x = round(x2); landmarks(5)._y = round(y2);
		landmarks(8)._x = round(x1); landmarks(8)._y = round(y1);
		//landmarks(1)._x = round(x2); landmarks(1)._y = round(y2);
		//landmarks(3)._x = round(x1); landmarks(3)._y = round(y1);
	}else{
		landmarks(5)._x = round(x1); landmarks(5)._y = round(y1);
		landmarks(8)._x = round(x2); landmarks(8)._y = round(y2);
		//landmarks(1)._x = round(x1); landmarks(1)._y = round(y1);
		//landmarks(3)._x = round(x2); landmarks(3)._y = round(y2);
	}
    //regenerate 3 1
	K = (landmarks(6)._y - landmarks(8)._y)/(landmarks(6)._x - landmarks(8)._x);
	K = K + M_PI/4;
    solveequation( x1,  x2,  y1,  y2,  Cx,  R,  K,  Cy);

	kp = (Cx-landmarks(6)._x)*(y1-landmarks(6)._y) - (Cy-landmarks(6)._y)*(x1-landmarks(6)._x);
	if(kp > 0){
		landmarks(1)._x = round(x2); landmarks(1)._y = round(y2);
		landmarks(3)._x = round(x1); landmarks(3)._y = round(y1);
	}else{
		landmarks(1)._x = round(x1); landmarks(1)._y = round(y1);
		landmarks(3)._x = round(x2); landmarks(3)._y = round(y2);
	}

	//regenerate 7 10
	K = -(landmarks(6)._x - landmarks(8)._x)/(landmarks(6)._y - landmarks(8)._y);
    solveequation( x1,  x2,  y1,  y2,  Cx,  R,  K,  Cy);

	kp = (Cx-landmarks(6)._x)*(y1-landmarks(6)._y) - (Cy-landmarks(6)._y)*(x1-landmarks(6)._x);
	if(kp > 0){
		//landmarks(1)._x = round(x2); landmarks(1)._y = round(y2);
		//landmarks(3)._x = round(x1); landmarks(3)._y = round(y1);
		landmarks(10)._x = round(x2); landmarks(10)._y = round(y2);
		landmarks(7)._x = round(x1); landmarks(7)._y = round(y1);
	}else{
        //landmarks(1)._x = round(x1); landmarks(1)._y = round(y1);
		//landmarks(3)._x = round(x2); landmarks(3)._y = round(y2);
		landmarks(10)._x = round(x1); landmarks(10)._y = round(y1);
		landmarks(7)._x = round(x2); landmarks(7)._y = round(y2);
	}
	//regenerate 2 4
	K = -(landmarks(6)._x - landmarks(8)._x)/(landmarks(6)._y - landmarks(8)._y);
	K = K + M_PI/4;
    solveequation( x1,  x2,  y1,  y2,  Cx,  R,  K,  Cy);
	
	kp = (Cx-landmarks(6)._x)*(y1-landmarks(6)._y) - (Cy-landmarks(6)._y)*(x1-landmarks(6)._x);
	if(kp > 0){
		landmarks(4)._x = round(x2); landmarks(4)._y = round(y2);
		landmarks(2)._x = round(x1); landmarks(2)._y = round(y1);
	}else{
        landmarks(4)._x = round(x1); landmarks(4)._y = round(y1);
		landmarks(2)._x = round(x2); landmarks(2)._y = round(y2);
	}
	for( i=0; i<landmarks.Size(); i++){
		target.ImageToWorld(landmarks(i));
	}
	landmarks.WriteVTK("surface\\newlandmarks.vtk");
	for( i=0; i<landmarks.Size(); i++){
		target.WorldToImage(landmarks(i));
	}
}

void irtkBep::Finalize(){
	
}

void irtkBep::EvaluateWallaverage(){
	int i,j,k;
	double integral = 0;
	int count = 0;
	for ( k = 0; k < compare.GetZ(); k++){
			for ( j = 0; j< compare.GetY(); j++){
				for ( i = 0; i<compare.GetX(); i++){
					if(threshold.GetAsDouble(i,j,k) == 3){
						integral += compare.GetAsDouble(i,j,k);
						count = count + 1;
					}
				}
			}
		}
}

void irtkBep::EvaluateWallThreshold(){
	int i,j,k;
	for ( k = 0; k < compare.GetZ(); k++){
			for ( j = 0; j< compare.GetY(); j++){
				for ( i = 0; i<compare.GetX(); i++){
					if(segmentation.GetAsDouble(i,j,k) > 0){
					}else
						 compare.PutAsDouble(i,j,k,0);
				}
			}
		}
	wallthreshold.Initialize(compare.GetImageAttributes());
	cout << "Blurring for segmentation ... "; cout.flush();
	irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(2, 0);
	blurring.SetInput (&compare);
	blurring.SetOutput(&compare);
	blurring.Run();
	cout << "done" << endl;
	int n = 2;
	double *mean, *var, *c;
	irtkRealPixel min, max;
	// Default parameters
	double treshold = 0.0001;
	int iterations = 50;
	//Default settings for Gaussian Mixture parameters
	mean = new double[n];
	var  = new double[n];
	c    = new double[n];

	compare.GetMinMax(&min,&max);

	for (i=0;i<n;i++)  mean[i] = min + i*(max-min)/(double) n;
	for (i=0;i<n;i++)  var[i] = ((max-min)/(double) n)*((max-min)/(double) n);
	for (i=0;i<n;i++)  c[i] = 1.0/(double) n;

	irtkEMClassification classification;
	classification.SetInput(compare);
	classification.SetPadding(0);
	classification.CreateMask();
	classification.InitialiseGMMParameters(n,mean,var,c);

	double rel_diff;
	i=1;
	cout << "Segmenting tissues: type Background Tissue Blood"<<endl;
	do {
		///cout << "Iteration = " << i << " / " << iterations << endl;
		rel_diff = classification.IterateGMM(i);
		i++;
	} while ((rel_diff>treshold)&&(i<iterations));
	cout << "Segmentation done, threshold matrix outputted"<<endl;
	
	irtkRealImage test1;
	classification.GetProbMap(1,test1);
	irtkRealImage test0;
	classification.GetProbMap(0,test0);

	/*for ( k = 0; k < test1.GetZ(); k++){
			for ( j = 0; j< test1.GetY(); j++){
				for ( i = 0; i<test1.GetX(); i++){
					if(compare.GetAsDouble(i,j,k) <= 0)
						 test1.PutAsDouble(i,j,k,0);
				}
			}
		}

	irtkImageGraphCut<float> graphcut;
	graphcut.SetInput(&compare,&test0,&test1);
	graphcut.SetOutput(&wallthreshold);
	graphcut.SetMode(2);
	graphcut.Run();*/

	classification.ConstructSegmentationNoBG(wallthreshold);
	for ( k = 0; k < compare.GetZ(); k++){
			for ( j = 0; j< compare.GetY(); j++){
				for ( i = 0; i<compare.GetX(); i++){
					if(compare.GetAsDouble(i,j,k) >= 
						(mean[1] - ((mean[1] - mean[0])*var[1])/
						((var[1] + var[0])*5))){
					}else if(wallthreshold.GetAsDouble(i,j,k)==2){
						 wallthreshold.PutAsDouble(i,j,k,1);
					}
				}
			}
		}

	delete mean;
	delete var;
	delete c;	
	char buffer[255];
	sprintf(buffer, "%swthres.nii",this->outputfilename);
	this->wallthreshold.Write(buffer);
}

void irtkBep::Bullseyeplot(int average){
  int i,t,swapped = 0;
  bep=new double*[compare.GetT()];                                                          
  for(i=0;i<compare.GetT();i++)             
	  bep[i]=new double[17];

  if(average == -1){
	  EvaluateWallThreshold();
	  EvaluateWallaverage();
  }

  for(t=0;t<compare.GetT();t++){
	  cout << "Evaluating " << t <<"th frame..";
	  
	  //Evaluate Apex
	  {
		  /*Sector apex;
		  apex.cx = landmarks(0)._x; apex.cy = landmarks(0)._y + compare.GetY()/4;
		  apex.x1 = landmarks(0)._x - compare.GetX()/4; apex.y1 = landmarks(0)._y - compare.GetY()/4;
		  apex.x2 = landmarks(0)._x + compare.GetX()/4; apex.y2 = landmarks(0)._y - compare.GetY()/4;
		  if(tip > aa)
			bep[t][16] = stintegral(compare, apex, aa, tip, t, BepMode);
		  else
			bep[t][16] = stintegral(compare, apex, tip, aa, t, BepMode);
			*/
		  bep[t][16] = 0;
	  }
		

	  //Evaluate Apical
	  {
		  Sector apical;
		  apical.cx = landmarks(3)._x;
		  apical.cy = landmarks(3)._y;
		  if ( am >  aa){		  swap(aa,am); swapped = 1;}
		  bep[t][12] = evaluatesector(compare, apical, landmarks(1), landmarks(2), am, aa,12,t,average);
		  bep[t][13] = evaluatesector(compare, apical, landmarks(2), landmarks(3), am, aa,13,t,average);
		  bep[t][14] = evaluatesector(compare, apical, landmarks(3), landmarks(4), am, aa,14,t,average);
		  bep[t][15] = evaluatesector(compare, apical, landmarks(4), landmarks(1), am, aa,15,t,average);
		  if ( swapped == 1){		  swap(aa,am);swapped = 0;}
	  }

	  //Evaluate Mid_Cavity
	  {
		  Sector mid;
		  mid.cx = landmarks(0)._x;
		  mid.cy = landmarks(0)._y;
		  if ( mb >  am){		  swap(mb,am); swapped = 1;}
		  bep[t][6] = evaluatesector(compare, mid, landmarks(5), landmarks(6), mb, am,6,t,average);
		  bep[t][7] = evaluatesector(compare, mid, landmarks(6), landmarks(7), mb, am,7,t,average);
		  bep[t][8] = evaluatesector(compare, mid, landmarks(7), landmarks(8), mb, am,8,t,average);
		  bep[t][9] = evaluatesector(compare, mid, landmarks(8), landmarks(9), mb, am,9,t,average);
		  bep[t][10] = evaluatesector(compare, mid, landmarks(9), landmarks(10), mb, am,10,t,average);
		  bep[t][11] = evaluatesector(compare, mid, landmarks(10), landmarks(5), mb, am,11,t,average);
		  if ( swapped == 1){		  swap(mb,am);swapped = 0;}		 
	  }

	  //Evaluate Basal
	  {
		  Sector basal;
		  basal.cx = landmarks(4)._x;
		  basal.cy = landmarks(4)._y;
		  if ( bottom >  mb){		  swap(mb,bottom); swapped = 1;}
		  bep[t][0] = evaluatesector(compare, basal, landmarks(5), landmarks(6), bottom, mb,0,t,average);
		  bep[t][1] = evaluatesector(compare, basal, landmarks(6), landmarks(7), bottom, mb,1,t,average);
		  bep[t][2] = evaluatesector(compare, basal, landmarks(7), landmarks(8), bottom, mb,2,t,average);
		  bep[t][3] = evaluatesector(compare, basal, landmarks(8), landmarks(9), bottom, mb,3,t,average);
		  bep[t][4] = evaluatesector(compare, basal, landmarks(9), landmarks(10), bottom, mb,4,t,average);
		  bep[t][5] = evaluatesector(compare, basal, landmarks(10), landmarks(5), bottom, mb,5,t,average);
		  if ( swapped == 1){		  swap(mb,bottom);swapped = 0;}		  
	  }
  }
   //Regulate the bulleye plot
	  double max = 0;
	  for(t=0;t<compare.GetT()-1;t++){
		  for(i=0; i<17; i++){
			  if(max < abs(bep[t][i])){
				  max = abs(bep[t][i]);
			  }
		  }
	  }	
	  //Output max to file
	  ofstream maxout(maxfilename,ios::app);
	  maxout<<max<<" ";

	  max = 0;
	  for(t=compare.GetT()-1;t<compare.GetT();t++){
		  for(i=0; i<17; i++){
			  if(max < abs(bep[t][i])){
				  max = abs(bep[t][i]);
			  }
		  }
	  }
	  maxout<<max<<endl;
	  maxout.close();
	 


	  //Output bulleyeplot to file
	  ofstream fout(outputfilename,ios::app);
	  for(t=0;t<compare.GetT();t++){
		  for(i=0; i<17; i++)
			  fout << bep[t][i] <<" ";
		  fout << endl;
	  }
	  fout.close();

	  for(i=0;i<compare.GetT();i++)
	  {
		  delete []bep[i];
		  bep[i]=NULL;
	  }
	  delete []bep;
	  bep=NULL;
	  cout<<endl;

}

void irtkBep::GenerateSegmentation(char* surfacefilename){
  int i,swapped = 0;
  irtkImageAttributes atr;
  atr = target.GetImageAttributes();
  atr._t = 1;
  segmentation.Initialize(atr);
  irtkRealImage radius(atr);
  
  	  
  //Evaluate Apex
  /*{
	  Sector apex;
	  apex.cx = landmarks(0)._x; apex.cy = landmarks(0)._y + target.GetY()/4;
	  apex.x1 = landmarks(0)._x - target.GetX()/4; apex.y1 = landmarks(0)._y - target.GetY()/4;
	  apex.x2 = landmarks(0)._x + target.GetX()/4; apex.y2 = landmarks(0)._y - target.GetY()/4;
	  if ( aa >  tip){		  swap(aa,tip); swapped = 1;}
	  for(k = aa; k < tip; k++)
		  for(j = 0; j<segmentation.GetY(); j++)
			  for(i=0;i<segmentation.GetX();i++){
				  if(isinside(i,j,apex) && threshold.GetAsDouble(i,j,k) == 3){
					  segmentation.PutAsDouble(i,j,k,17);
				  }
			  }
	 if ( swapped == 1){		  swap(aa,tip);swapped = 0;}
  }*/


	  //Evaluate Apical
	  {
		  Sector apical;
		  apical.cx = olandmarks(3)._x;
		  apical.cy = olandmarks(3)._y;
		  if ( am >  aa){		  swap(aa,am); swapped = 1;}
		  segmentsector(radius,segmentation, apical, landmarks(1), landmarks(2), am, aa,13);
		  segmentsector(radius,segmentation, apical, landmarks(2), landmarks(3), am, aa,14);
		  segmentsector(radius,segmentation, apical, landmarks(3), landmarks(4), am, aa,15);
		  segmentsector(radius,segmentation, apical, landmarks(4), landmarks(1), am, aa,16);
		  if ( swapped == 1){		  swap(aa,am);swapped = 0;}
	  }

	  //Evaluate Mid_Cavity
	  {
		  Sector mid;
		  mid.cx = olandmarks(0)._x;
		  mid.cy = olandmarks(0)._y;
		  if ( mb >  am){		  swap(mb,am); swapped = 1;}
		  segmentsector(radius,segmentation, mid, landmarks(5), landmarks(6), mb, am,7);
		  segmentsector(radius,segmentation, mid, landmarks(6), landmarks(7), mb, am,8);
		  segmentsector(radius,segmentation, mid, landmarks(7), landmarks(8), mb, am,9);
		  segmentsector(radius,segmentation, mid, landmarks(8), landmarks(9), mb, am,10);
		  segmentsector(radius,segmentation, mid, landmarks(9), landmarks(10), mb, am,11);
		  segmentsector(radius,segmentation, mid, landmarks(10), landmarks(5), mb, am,12);
		  if ( swapped == 1){		  swap(mb,am);swapped = 0;}
	  }

	  //Evaluate Basal
	  {
		  Sector basal;
		  basal.cx = olandmarks(4)._x;
		  basal.cy = olandmarks(4)._y;
		  if ( bottom >  mb){		  swap(mb,bottom); swapped = 1;}
		  segmentsector(radius,segmentation, basal, landmarks(5), landmarks(6), bottom, mb,1);
		  segmentsector(radius,segmentation, basal, landmarks(6), landmarks(7), bottom, mb,2);
		  segmentsector(radius,segmentation, basal, landmarks(7), landmarks(8), bottom, mb,3);
		  segmentsector(radius,segmentation, basal, landmarks(8), landmarks(9), bottom, mb,4);
		  segmentsector(radius,segmentation, basal, landmarks(9), landmarks(10), bottom, mb,5);
		  segmentsector(radius,segmentation, basal, landmarks(10), landmarks(5), bottom, mb,6);
		  if ( swapped == 1){		  swap(mb,bottom);swapped = 0;}
	  }
	  char buffer[255];
	  sprintf(buffer, "%s\\segmentation.nii",surfacefilename);	  
      segmentation.Write(buffer);
	  sprintf(buffer, "%s\\radius.nii",surfacefilename);
      radius.Write(buffer);
	  sprintf(buffer, "%s\\centerpoints.vtk",surfacefilename);
	  for( i=0; i<centerpoints.Size(); i++){
		segmentation.ImageToWorld(centerpoints(i));
	  }
	  centerpoints.WriteVTK(buffer);
}
	
void irtkBep::Compareoverall(){
	  irtkGreyImage tinte,sinte;
      int i,j,k,t;
	  irtkImageAttributes atr;
	  //find out number of time frame
	  atr = target.GetImageAttributes();

	  compare.Initialize(atr);
	  cout << "Evaluating original data "<< endl;
	  analysis(osequence,tinte);
	  cout << "Done..." << endl;
	  cout << "Evaluating 1yr later data "<< endl;
	  analysis(lsequence,sinte);
	  cout << "Done..." << endl;

	  //correspond sinte to tinte

	  for(t =0; t < atr._t; t++){
		  for(k = 0; k < atr._z; k++)
			  for(j = 0; j < atr._y; j++)
				  for(i = 0; i < atr._x; i++){
					  double tmp = (sinte.GetAsDouble(i,j,k,t) - tinte.GetAsDouble(i,j,k,t));
					  if(tmp != 0){
						  //double norm = abs(tinte.GetAsDouble(i,j,k,t) + sinte.GetAsDouble(i,j,k,t));
						  //tmp = tmp/norm;
						  compare.PutAsDouble(i,j,k,t,tmp);	
					  }else{
						  compare.PutAsDouble(i,j,k,t,0);
					  }
				  }
	  }

	  compare.Write("surface\\idcompare.nii");
	  //compare.Write("surface\\idcompare.vtk");
	 
	 // tinte.Write("tinte.nii");
	 // sinte.Write("sinte.nii");
    /*cout << "Blurring compare ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(1, 1);
    blurring.SetInput (&compare);
    blurring.SetOutput(&compare);
    blurring.Run();
    cout << "done" << endl;*/

}

void irtkBep::Comparedifference(){
	irtkRealImage inte;
      int i,j,k,t;
	  irtkImageAttributes atr;
	  //find out number of time frame
	  atr = target.GetImageAttributes();
	  compare.Initialize(atr);
	  cout << "Evaluating data "<< endl;
	  analysis(inte);
	  cout << "Done..." << endl;

	  for(t =1; t < atr._t; t++){
		  for(k = 0; k < atr._z; k++)
			  for(j = 0; j < atr._y; j++)
				  for(i = 0; i < atr._x; i++){
					  compare.PutAsDouble(i,j,k,t-1,inte.GetAsDouble(i,j,k,t));
					  compare.PutAsDouble(i,j,k,atr._t-1,compare.GetAsDouble(i,j,k,atr._t-1)+inte.GetAsDouble(i,j,k,t));
				  }
	  }
	 
    /*cout << "Blurring compare ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(1, 1);
    blurring.SetInput (&compare);
    blurring.SetOutput(&compare);
    blurring.Run();
    cout << "done" << endl;*/
    compare.Write("surface\\dmcompare.nii");
	//compare.Write("surface\\dmcompare.vtk");

}

void irtkBep::Compareradius(){
	  irtkRealImage inte;
      int i,j,k,t;
	  irtkImageAttributes atr;
	  //find out number of time frame
	  atr = target.GetImageAttributes();

	  compare.Initialize(atr);
	  cout << "Evaluating data "<< endl;
	  analysis(inte,landmarks(0));
	  cout << "Done..." << endl;

	  for(t =1; t < atr._t; t++){
		  for(k = 0; k < atr._z; k++)
			  for(j = 0; j < atr._y; j++)
				  for(i = 0; i < atr._x; i++){
					  compare.PutAsDouble(i,j,k,t-1,inte.GetAsDouble(i,j,k,t));
					  compare.PutAsDouble(i,j,k,atr._t-1,compare.GetAsDouble(i,j,k,atr._t-1)+inte.GetAsDouble(i,j,k,t));
				  }
	  }

	  /*cout << "Blurring compare ... "; cout.flush();
	  irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(1, 1);
	  blurring.SetInput (&compare);
	  blurring.SetOutput(&compare);
	  blurring.Run();
	  cout << "done" << endl;*/
	  compare.Write("surface\\rdcompare.nii");
	  //compare.Write("surface\\rdcompare.vtk");
}

void irtkBep::analysis(irtkRealImage& integral, irtkPoint& center){
	double i,j,k,x,y,z;
	double x1,y1,z1;
	double d1,d2;
	int l,t;
	irtkRealImage sub;
	irtkImageAttributes atr;
	// Create initial multi-level free-form deformation
	irtkMultiLevelFreeFormTransformation *mffd1,*mffd2;
	//find out number of time frame
	atr = target.GetImageAttributes();
	t = atr._t;
	atr._t = 1;

	//generate integral image
	sub.Initialize(atr);

	//compute original radius using sub
	for(k = 0; k < atr._z; k++)
		for(j = 0; j < atr._y; j++)
			for(i = 0; i < atr._x; i++){
				if( round(threshold.GetAsDouble(i,j,k)) == 3){
					x = center._x; y = center._y;
					x = i - x; y = j - y;
					x = x * sub.GetXSize(); y = y * sub.GetYSize();
					sub.PutAsDouble(i,j,k,sqrt(x*x+y*y));
				}
			}
			//generate integral
			atr._t = t;
			integral.Initialize(atr);

			//loop over the target image
			for(l = 1; l < t; l++){
				mffd1 = NULL; mffd2 = NULL;
				//read in transformation
				char buffer[255];
				sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", osequence, 0, l);

				irtkTransformation *transform = irtkTransformation::New(buffer);
				if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
					mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
				} else {
					if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
						mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
					} else {
						if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
							mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
						} else {
							cerr << "Input transformation is not of type rigid or affine " << endl;
							cerr << "or multi-level free form deformation" << endl;
							exit(1);
				  }
			  }
				}
				delete transform;

				//read another transformation
				sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", lsequence, 0, l);

				transform = irtkTransformation::New(buffer);
				if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
					mffd2 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
				} else {
					if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
						mffd2 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
					} else {
						if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
							mffd2 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
						} else {
							cerr << "Input transformation is not of type rigid or affine " << endl;
							cerr << "or multi-level free form deformation" << endl;
							exit(1);
				  }
			  }
				}
				delete transform;

				//loop over voxels
				for(k = 0; k < atr._z; k++)
					for(j = 0; j < atr._y; j++)
						for(i = 0; i < atr._x; i++){
							if( round(threshold.GetAsDouble(i,j,k)) == 3){
								//evaluate displacement add to integral image
								//if(i==42&&j==53&&k==4){
								//  int text = 1;
						  //}
								x = i; y = j; z = k;
								x1 = i; y1 = j; z1 = k;
								sub.ImageToWorld(x,y,z);
								sub.ImageToWorld(x1,y1,z1);
								mffd1->Transform(x,y,z);
								mffd2->Transform(x1,y1,z1);
								sub.WorldToImage(x,y,z);
								sub.WorldToImage(x1,y1,z1);
								x = x - center._x; y = y - center._y;
								x1 = x1 - center._x; y1 = y1 - center._y;
								x = x * sub.GetXSize(); y = y * sub.GetYSize(); z = z * sub.GetZSize();
								x1 = x1 * sub.GetXSize(); y1 = y1 * sub.GetYSize(); z1 = z1 * sub.GetZSize();
								d1 = sub.GetAsDouble(i,j,k) - sqrt(x*x + y*y);
								d2 = sub.GetAsDouble(i,j,k) - sqrt(x1*x1 + y1*y1);
								double tmp;
								if ((d1+d2) != 0){
									tmp = (d2-d1);
								}else{
									tmp = 0;
						  }
								integral.PutAsDouble(i,j,k,l,tmp);
					  }
				  }

						delete mffd1;
						delete mffd2;
			}
}

void irtkBep::AnalysisRadius(int prefix){
	double i,j,k,x,y,z;
	double d1;
	int l,t,swapped;
	swapped = 0;
	irtkRealImage sub;
	irtkImageAttributes atr;
	// Create initial multi-level free-form deformation
	irtkMultiLevelFreeFormTransformation *mffd1;
	//find out number of time frame
	atr = target.GetImageAttributes();
	t = atr._t;
	atr._t = 1;

	//generate integral image
	sub.Initialize(atr);
	//evaluate initial radial distance
	//apical
	if ( am >  aa){		  swap(aa,am); swapped = 1;}
	for(k = am; k < aa; k++){
		for(j = 0; j < atr._y; j++){
			for(i = 0; i < atr._x; i++){
				if( round(threshold.GetAsDouble(i,j,k)) == 3){
					x = landmarks(3)._x; y = landmarks(3)._y;
					x = i - x; y = j - y;
					x = x * sub.GetXSize(); y = y * sub.GetYSize();
					sub.PutAsDouble(i,j,k,sqrt(x*x+y*y));
				}
			}
		}
	}
	if ( swapped == 1){		  swap(aa,am);swapped = 0;}
	//middle
	if ( mb >  am){		  swap(mb,am); swapped = 1;}
	for(k = mb; k < am; k++){
		for(j = 0; j < atr._y; j++){
			for(i = 0; i < atr._x; i++){
				if( round(threshold.GetAsDouble(i,j,k)) == 3){
					x = landmarks(0)._x; y = landmarks(0)._y;
					x = i - x; y = j - y;
					x = x * sub.GetXSize(); y = y * sub.GetYSize();
					sub.PutAsDouble(i,j,k,sqrt(x*x+y*y));
				}
			}
		}
	}
	if ( swapped == 1){		  swap(mb,am);swapped = 0;}
	//base
	if ( bottom >  mb){		  swap(mb,bottom); swapped = 1;}
	for(k = bottom; k < mb; k++){
		for(j = 0; j < atr._y; j++){
			for(i = 0; i < atr._x; i++){
				if( round(threshold.GetAsDouble(i,j,k)) == 3){
					x = landmarks(4)._x; y = landmarks(4)._y;
					x = i - x; y = j - y;
					x = x * sub.GetXSize(); y = y * sub.GetYSize();
					sub.PutAsDouble(i,j,k,sqrt(x*x+y*y));
				}
			}
		}
	}
	if ( swapped == 1){		  swap(mb,bottom);swapped = 0;}
	//generate integral
	atr._t = t;
	compare.Initialize(atr);

	//loop over the target image
	for(l = 1; l < t; l++){
		mffd1 = NULL;
		//read in transformation
		char buffer[255];
		sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", osequence, prefix, l);

		irtkTransformation *transform = irtkTransformation::New(buffer);
		if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
			mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
		} else {
			if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
				mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
			} else {
				if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
					mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
				} else {
					cerr << "Input transformation is not of type rigid or affine " << endl;
					cerr << "or multi-level free form deformation" << endl;
					exit(1);
				}
			}
		}
		delete transform;
		
		//evaluate initial radial distance
		//apical
		if ( am >  aa){		  swap(aa,am); swapped = 1;}
		for(k = am; k < aa; k++){
			for(j = 0; j < atr._y; j++){
				for(i = 0; i < atr._x; i++){
					if( round(threshold.GetAsDouble(i,j,k)) == 3){
						x = i; y = j; z = k;
						sub.ImageToWorld(x,y,z);
						mffd1->Transform(x,y,z);
						sub.WorldToImage(x,y,z);
						x = x - landmarks(3)._x; y = y - landmarks(3)._y;
						x = x * sub.GetXSize(); y = y * sub.GetYSize();
						d1 = sub.GetAsDouble(i,j,k) - sqrt(x*x + y*y);
						compare.PutAsDouble(i,j,k,l-1,d1);
						compare.PutAsDouble(i,j,k,t-1,d1+compare.GetAsDouble(i,j,k,t-1));
					}
				}
			}
		}
		if ( swapped == 1){		  swap(aa,am);swapped = 0;}
		//middle
		if ( mb >  am){		  swap(mb,am); swapped = 1;}
		for(k = mb; k < am; k++){
			for(j = 0; j < atr._y; j++){
				for(i = 0; i < atr._x; i++){
					if( round(threshold.GetAsDouble(i,j,k)) == 3){
						x = i; y = j; z = k;
						sub.ImageToWorld(x,y,z);
						mffd1->Transform(x,y,z);
						sub.WorldToImage(x,y,z);
						x = x - landmarks(0)._x; y = y - landmarks(0)._y;
						x = x * sub.GetXSize(); y = y * sub.GetYSize();
						d1 = sub.GetAsDouble(i,j,k) - sqrt(x*x + y*y);
						compare.PutAsDouble(i,j,k,l-1,d1);
						compare.PutAsDouble(i,j,k,t-1,d1+compare.GetAsDouble(i,j,k,t-1));
					}
				}
			}
		}
		if ( swapped == 1){		  swap(mb,am);swapped = 0;}
		//base
		if ( bottom >  mb){		  swap(mb,bottom); swapped = 1;}
		for(k = bottom; k < mb; k++){
			for(j = 0; j < atr._y; j++){
				for(i = 0; i < atr._x; i++){
					if( round(threshold.GetAsDouble(i,j,k)) == 3){
						x = i; y = j; z = k;
						sub.ImageToWorld(x,y,z);
						mffd1->Transform(x,y,z);
						sub.WorldToImage(x,y,z);
						x = x - landmarks(4)._x; y = y - landmarks(4)._y;
						x = x * sub.GetXSize(); y = y * sub.GetYSize();
						d1 = sub.GetAsDouble(i,j,k) - sqrt(x*x + y*y);
						compare.PutAsDouble(i,j,k,l-1,d1);
						compare.PutAsDouble(i,j,k,t-1,d1+compare.GetAsDouble(i,j,k,t-1));
					}
				}
			}
		}
		if ( swapped == 1){		  swap(mb,bottom);swapped = 0;}

		delete mffd1;
	}
	compare.Write("surface\\rdmotion.nii");
}

void irtkBep::AnalysisMotion(int prefix){
	double i,j,k,x,y,z;
	double x1,y1,z1;
	double d1;
	int l,t,swapped;
	swapped = 0;
	irtkImageAttributes atr;
	// Create initial multi-level free-form deformation
	irtkMultiLevelFreeFormTransformation *mffd1;
	//find out number of time frame
	atr = target.GetImageAttributes();
	t = atr._t;
	compare.Initialize(atr);

	//loop over the target image
	for(l = 1; l < t; l++){
		mffd1 = NULL;
		//read in transformation
		char buffer[255];
		sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", osequence, prefix, l);

		irtkTransformation *transform = irtkTransformation::New(buffer);
		if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
			mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
		} else {
			if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
				mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
			} else {
				if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
					mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
				} else {
					cerr << "Input transformation is not of type rigid or affine " << endl;
					cerr << "or multi-level free form deformation" << endl;
					exit(1);
				}
			}
		}
		delete transform;
		
		//evaluate initial displacement
		for(k = 0; k < atr._z; k++){
			for(j = 0; j < atr._y; j++){
				for(i = 0; i < atr._x; i++){
					if( round(segmentation.GetAsDouble(i,j,k)) > 0){
						x = i; y = j; z = k;
						segmentation.ImageToWorld(x,y,z);
						x1 = x; y1 = y; z1 = z;
						mffd1->Transform(x,y,z);
						x = x - x1; y = y - y1; z = z - z1;
						d1 = sqrt(x*x + y*y + z*z);
						compare.PutAsDouble(i,j,k,l-1,d1);
						compare.PutAsDouble(i,j,k,t-1,d1+compare.GetAsDouble(i,j,k,t-1));
					}
				}
			}
		}
		delete mffd1;
	}
	compare.Write("surface\\motion.nii");
}

void irtkBep::AnalysisStrain(int prefix,int mode){
	double i,j,k,x,y,z,d1;
	double tx[3],ty[3],tz[3];
	double v[3],u[3],w[3],tmp[3];
	int m,l,t,swapped;
	swapped = 0;
	d1 = 0;
	irtkImageAttributes atr;
	// Create initial multi-level free-form deformation
	irtkMultiLevelFreeFormTransformation *mffd1;
	//find out number of time frame
	atr = target.GetImageAttributes();
	t = atr._t;
	compare.Initialize(atr);

	//loop over the target image
	for(l = 1; l < t; l++){
		mffd1 = NULL;
		//read in transformation
		char buffer[255];
		sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", osequence, prefix, l);

		irtkTransformation *transform = irtkTransformation::New(buffer);
		if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
			mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
		} else {
			if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
				mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
			} else {
				if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
					mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
				} else {
					cerr << "Input transformation is not of type rigid or affine " << endl;
					cerr << "or multi-level free form deformation" << endl;
					exit(1);
				}
			}
		}
		delete transform;
		
		//evaluate initial strain
		for(k = 1; k < atr._z-1; k++){
			for(j = 1; j < atr._y-1; j++){
				for(i = 1; i < atr._x-1; i++){
					if( round(segmentation.GetAsDouble(i,j,k)) > 0){
						if( mode == 0){
							x = i-0.5; y = j-0.5; z = k-0.5;
							tx[0] = x+1; ty[0] = y; tz[0] = z;
							tx[1] = x; ty[1] = y+1; tz[1] = z;
							tx[2] = x; ty[2] = y; tz[2] = z+1;
							segmentation.ImageToWorld(x,y,z);						
							for(m = 0; m < 3; m++){
								segmentation.ImageToWorld(tx[m],ty[m],tz[m]);							
							}
							//uyvz - uzvy, uzvx - uxvz, uxvy - uyvx
							//uxux+uyuy+uzuz
							u[0] = tx[0] - x; u[1] = ty[0] - y; u[2] = tz[0] - z;
							v[0] = tx[1] - x; v[1] = ty[1] - y; v[2] = tz[1] - z;
							w[0] = tx[2] - x; w[1] = ty[2] - y; w[2] = tz[2] - z;
							tmp[0] = u[1]*v[2] - u[2]*v[1];
							tmp[1] = u[2]*v[0] - u[0]*v[2];
							tmp[2] = u[0]*v[1] - u[1]*v[0];
							d1 = fabs(w[0]*tmp[0] + w[1]*tmp[1] + w[2]*tmp[2]);
							mffd1->Transform(x,y,z);
							for(m = 0; m < 3; m++){
								mffd1->Transform(tx[m],ty[m],tz[m]);
							}
							u[0] = tx[0] - x; u[1] = ty[0] - y; u[2] = tz[0] - z;
							v[0] = tx[1] - x; v[1] = ty[1] - y; v[2] = tz[1] - z;
							w[0] = tx[2] - x; w[1] = ty[2] - y; w[2] = tz[2] - z;
							tmp[0] = u[1]*v[2] - u[2]*v[1];
							tmp[1] = u[2]*v[0] - u[0]*v[2];
							tmp[2] = u[0]*v[1] - u[1]*v[0];
							d1 = fabs(w[0]*tmp[0] + w[1]*tmp[1] + w[2]*tmp[2])/d1;
						}else if(mode == 1){
							// longitudinal strain
							x = i-0.5; y = j-0.5; z = k-0.5;
							tx[0] = x+1; ty[0] = y; tz[0] = z;
							tx[1] = x; ty[1] = y+1; tz[1] = z;
							tx[2] = x; ty[2] = y; tz[2] = z+1;
							segmentation.ImageToWorld(x,y,z);						
							for(m = 0; m < 3; m++){
								segmentation.ImageToWorld(tx[m],ty[m],tz[m]);							
							}
							mffd1->Transform(x,y,z);
							segmentation.WorldToImage(x,y,z);
							for(m = 0; m < 3; m++){
								mffd1->Transform(tx[m],ty[m],tz[m]);
								segmentation.WorldToImage(tx[m],ty[m],tz[m]);
							}
							//w[0] = tx[2] - x; w[1] = ty[2] - y; 
							w[2] = tz[2] - z;
							d1 = sqrt(//w[0]*w[0] + w[1]*w[1] + 
								w[2]*w[2]);
						}else if(mode == 2){
							// radial and circumferential strain
							x = i-0.5; y = j-0.5; z = k-0.5;
							tx[0] = x+1; ty[0] = y; tz[0] = z;
							tx[1] = x; ty[1] = y+1; tz[1] = z;
							tx[2] = x; ty[2] = y; tz[2] = z+1;
							segmentation.ImageToWorld(x,y,z);						
							for(m = 0; m < 3; m++){
								segmentation.ImageToWorld(tx[m],ty[m],tz[m]);							
							}
							//uyvz - uzvy, uzvx - uxvz, uxvy - uyvx
							//uxux+uyuy+uzuz
							mffd1->Transform(x,y,z);
							segmentation.WorldToImage(x,y,z);
							for(m = 0; m < 3; m++){
								mffd1->Transform(tx[m],ty[m],tz[m]);
								segmentation.WorldToImage(tx[m],ty[m],tz[m]);
							}
							u[0] = tx[0] - x; u[1] = ty[0] - y; u[2] = tz[0] - z;
							v[0] = tx[1] - x; v[1] = ty[1] - y; v[2] = tz[1] - z;
							w[0] = tx[2] - x; w[1] = ty[2] - y; w[2] = tz[2] - z;
							//tmp[0] = u[1]*v[2] - u[2]*v[1];
							//tmp[1] = u[2]*v[0] - u[0]*v[2];
							tmp[2] = u[0]*v[1] - u[1]*v[0];
							d1 = sqrt(//tmp[0]*tmp[0] + tmp[1]*tmp[1] + 
								tmp[2]*tmp[2]);
						}else if(mode == 3){
							// reversed longitudinal global strain
							// radial and circumferential strain
							x = i-0.5; y = j-0.5; z = k-0.5;
							tx[0] = x+1; ty[0] = y; tz[0] = z;
							tx[1] = x; ty[1] = y+1; tz[1] = z;
							tx[2] = x; ty[2] = y; tz[2] = z+1;
							segmentation.ImageToWorld(x,y,z);						
							for(m = 0; m < 3; m++){
								segmentation.ImageToWorld(tx[m],ty[m],tz[m]);							
							}
							//uyvz - uzvy, uzvx - uxvz, uxvy - uyvx
							//uxux+uyuy+uzuz
							mffd1->Transform(x,y,z);
							segmentation.WorldToImage(x,y,z);
							for(m = 0; m < 3; m++){
								mffd1->Transform(tx[m],ty[m],tz[m]);
								segmentation.WorldToImage(tx[m],ty[m],tz[m]);
							}
							u[0] = tx[0] - x; u[1] = ty[0] - y; u[2] = tz[0] - z;
							v[0] = tx[1] - x; v[1] = ty[1] - y; v[2] = tz[1] - z;
							w[0] = tx[2] - x; w[1] = ty[2] - y; w[2] = tz[2] - z;
							//tmp[0] = u[1]*v[2] - u[2]*v[1];
							//tmp[1] = u[2]*v[0] - u[0]*v[2];
							tmp[2] = u[0]*v[1] - u[1]*v[0];
							d1 = sqrt(tmp[2]*tmp[2]);
							if(w[2] >0.01 || w[2] < -0.01){
								d1 = d1/sqrt(w[2]*w[2]);
								if(d1>2) d1 = 2;
							}else
								d1 = 2;
						}
						compare.PutAsDouble(i,j,k,l-1,d1);
						compare.PutAsDouble(i,j,k,t-1,d1+compare.GetAsDouble(i,j,k,t-1));
					}
				}
			}
		}
		delete mffd1;
	}
	compare.Write("surface\\strain.nii");
}

void irtkBep::analysis(irtkRealImage& integral){
	double i,j,k,x,y,z;
	double x1,y1,z1;
	double d1,d2;
	int l,t;
	irtkGreyImage sub;
	irtkImageAttributes atr;
	// Create initial multi-level free-form deformation
	irtkMultiLevelFreeFormTransformation *mffd1,*mffd2;
	//find out number of time frame
	atr = target.GetImageAttributes();
	t = atr._t;
	atr._t = 1;

	//generate integral image
	sub.Initialize(atr);
	atr._t = t;
	integral.Initialize(atr);

	//loop over the target image
	for(l = 1; l < t; l++){
		mffd1 = NULL; mffd2 = NULL;
		//read in transformation
		char buffer[255];
		sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", osequence, 0, l);

		irtkTransformation *transform = irtkTransformation::New(buffer);
		if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
			mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
		} else {
			if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
				mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
			} else {
				if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
					mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
				} else {
					cerr << "Input transformation is not of type rigid or affine " << endl;
					cerr << "or multi-level free form deformation" << endl;
					exit(1);
				}
			}
		}
		delete transform;

		//read another transformation
		sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", lsequence, 0, l);

		transform = irtkTransformation::New(buffer);
		if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
			mffd2 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
		} else {
			if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
				mffd2 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
			} else {
				if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
					mffd2 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
				} else {
					cerr << "Input transformation is not of type rigid or affine " << endl;
					cerr << "or multi-level free form deformation" << endl;
					exit(1);
				}
			}
		}
		delete transform;

		//loop over voxels
		for(k = 0; k < atr._z; k++)
			for(j = 0; j < atr._y; j++)
				for(i = 0; i < atr._x; i++){
					if( round(threshold.GetAsDouble(i,j,k)) == 3){
						//evaluate displacement add to integral image
						x = i; y = j; z = k;
						x1 = i; y1 = j; z1 = k;
						sub.ImageToWorld(x,y,z);
						sub.ImageToWorld(x1,y1,z1);
						mffd1->Transform(x,y,z);
						mffd2->Transform(x1,y1,z1);
						sub.WorldToImage(x,y,z);
						sub.WorldToImage(x1,y1,z1);
						x = x - i; y = y - j; z = z - k;
						x1 = x1 - i; y1 = y1 - j; z1 = z1 - k;
						x = x * sub.GetXSize(); y = y * sub.GetYSize(); z = z * sub.GetZSize();
						x1 = x1 * sub.GetXSize(); y1 = y1 * sub.GetYSize(); z1 = z1 * sub.GetZSize();
						d1 = sqrt(x*x + y*y + z*z);
						d2 = sqrt(x1*x1 + y1*y1 + z1*z1);
						double tmp;
						if ((d1+d2) != 0){
							tmp = (d2-d1);
						}else{
							tmp = 0;
						}
						integral.PutAsDouble(i,j,k,l,tmp);
					}
				}

				delete mffd1;
				delete mffd2;
	}
}

void irtkBep::analysis(char* sequence,irtkGreyImage& integral){
	double i,j,k,x,y,z;
	int l,t;
	irtkGreyImage sub;
	irtkImageAttributes atr;
	// Create initial multi-level free-form deformation
	irtkMultiLevelFreeFormTransformation *mffd;
	//find out number of time frame
	atr = target.GetImageAttributes();
	t = atr._t;
	atr._t = 1;

	//generate integral image
	sub.Initialize(atr);
	atr._t = t;
	integral.Initialize(atr);

	//loop over the target image
	for(l = 1; l < t; l++){
		mffd = NULL;
		//read in transformation
		char buffer[255];
		sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", sequence, 0, l);

		irtkTransformation *transform = irtkTransformation::New(buffer);
		if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
			mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
		} else {
			if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
				mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
			} else {
				if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
					mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
				} else {
					cerr << "Input transformation is not of type rigid or affine " << endl;
					cerr << "or multi-level free form deformation" << endl;
					exit(1);
				}
			}
		}
		delete transform;

		//loop over voxels
		for(k = 0; k < atr._z; k++)
			for(j = 0; j < atr._y; j++)
				for(i = 0; i < atr._x; i++){
					if( round(threshold.GetAsDouble(i,j,k)) == 3){
						//evaluate displacement add to integral image
						x = i; y = j; z = k;
						sub.ImageToWorld(x,y,z);
						mffd->Transform(x,y,z);
						sub.WorldToImage(x,y,z);
						x = x - i; y = y - j; z = z - k;
						x = x * sub.GetXSize(); y = y * sub.GetYSize(); z = z * sub.GetZSize();
						integral.PutAsDouble(i,j,k,l,integral.GetAsDouble(i,j,k,l-1) + sqrt(x*x + y*y + z*z));
					}
				}

				delete mffd;
	}
}

irtkBep::irtkBep(){
	outputfilename = NULL;
	osequence = NULL;
	lsequence = NULL;
}
