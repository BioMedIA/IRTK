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

#include <irtkRegistration2.h>

#include <irtkPatchBasedSegmentation.h>

char *cine_name = NULL;
char *seg_ED_name = NULL;
char *seg_ES_name = NULL;
char *output_name = NULL;
char *parin_name = NULL;

void usage()
{
	cerr << "Usage: cardiacpropergation [cine] [ED segmentation] [ES segmentation] [Segmentation output]" << endl;
	cerr << "<-parin file>        Read motion tracking parameter from file" << endl;
	cerr << "<-debug>             Debug mode" << endl;
	exit(1);
}

int main( int argc, char** argv )
{
	int ok,i,j,k,t,oldt,esphase,frames,debug;
	short cine_max,cine_min,cinedis;
	double *similarity,*smoothsimilarity,dif;
	// Check command line
	if (argc < 5) {
		usage();
	}

	// Parse source and target images
	cine_name = argv[1];
	argc--;
	argv++;
	seg_ED_name = argv[1];
	argc--;
	argv++;
	seg_ES_name = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;

	debug = false;

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-parin") == 0)) {
			argc--;
			argv++;
			ok = true;
			parin_name = argv[1];
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-debug") == 0)) {
			argc--;
			argv++;
			ok = true;
			debug = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	/// First identify ED and ES
	// Create images
	irtkGreyImage cine(cine_name);
	irtkGreyImage segmentation_ed(seg_ED_name);
	irtkGreyImage segmentation_es(seg_ES_name);
	irtkGreyImage output(cine.GetImageAttributes());
	irtkGreyImage propergated_segmentation_ed(segmentation_ed.GetImageAttributes());
	irtkGreyImage propergated_segmentation_es(segmentation_es.GetImageAttributes());
	irtkGreyImage tmp_segmentation(segmentation_ed.GetImageAttributes());
	irtkRealImage propergated_image_ed(segmentation_ed.GetImageAttributes());
	irtkRealImage propergated_image_es(segmentation_es.GetImageAttributes());
	irtkRealImage tmp_image(segmentation_ed.GetImageAttributes());
	irtkGreyImage blured;
	blured.Initialize(cine.GetImageAttributes());

	irtkGaussianBlurring<irtkGreyPixel> gaussianBlurring(1);
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

	/// Second, propergate from ED to all
	int index, x, y, z;
	/// Forward
	// create transformations and registrations
	irtkMultiLevelFreeFormTransformation **mffd_ed_f = new irtkMultiLevelFreeFormTransformation*[cine.GetT()];
	for(i = 0; i < cine.GetT(); i++){
		t = i;

		if(i > 0){
			mffd_ed_f[i] = new irtkMultiLevelFreeFormTransformation(*mffd_ed_f[i-1]);
		}else{
			mffd_ed_f[i] = new irtkMultiLevelFreeFormTransformation;
		}

		irtkImageFreeFormRegistration2 *registration = new irtkImageFreeFormRegistration2;

		// Combine images
		irtkImageAttributes attr = cine.GetImageAttributes();
		attr._t = 1;
		irtkGreyImage *target = new irtkGreyImage(attr);
		irtkGreyImage *source = new irtkGreyImage(attr);
		for (z = 0; z < target->GetZ(); z++) {
			for (y = 0; y < target->GetY(); y++) {
				for (x = 0; x < target->GetX(); x++) {
					target->Put(x, y, z, 0, cine.Get(x, y, z, 0));
					source->Put(x, y, z, 0, cine.Get(x, y, z, t));
				}
			}
		}

		// Set input and output for the registration filter
		registration->SetInput(target, source);
		registration->SetOutput(mffd_ed_f[i]);

		// Make an initial Guess for the parameters.
		registration->GuessParameter();
		// Overrride with any the user has set.
		if (parin_name != NULL) {
			registration->irtkImageRegistration2::Read(parin_name);
		}

		// Run registration filter
		registration->Run();

		if (registration->GetMFFDMode()) {
			mffd_ed_f[i]->CombineLocalTransformation();
		}

		if (debug) {
			char buffer[255];
			sprintf(buffer, "debug_ed_f_%.2d.dof.gz", i);
			mffd_ed_f[i]->irtkTransformation::Write(buffer);
		}

		delete registration;
		delete target;
		delete source;
	}
	/// Backward
	// create transformations and registrations
	irtkMultiLevelFreeFormTransformation **mffd_ed_b = new irtkMultiLevelFreeFormTransformation*[cine.GetT()];
	for(i = 0; i < cine.GetT(); i++){
		oldt = t;
		t = cine.GetT() - i;
		if(t == cine.GetT())
			t = 0;

		if(i > 0){
			mffd_ed_b[t] = new irtkMultiLevelFreeFormTransformation(*mffd_ed_b[oldt]);
		}else{
			mffd_ed_b[t] = new irtkMultiLevelFreeFormTransformation;
		}

		irtkImageFreeFormRegistration2 *registration = new irtkImageFreeFormRegistration2;

		// Combine images
		irtkImageAttributes attr = cine.GetImageAttributes();
		attr._t = 1;
		irtkGreyImage *target = new irtkGreyImage(attr);
		irtkGreyImage *source = new irtkGreyImage(attr);
		for (z = 0; z < target->GetZ(); z++) {
			for (y = 0; y < target->GetY(); y++) {
				for (x = 0; x < target->GetX(); x++) {
					target->Put(x, y, z, 0, cine.Get(x, y, z, 0));
					source->Put(x, y, z, 0, cine.Get(x, y, z, t));
				}
			}
		}

		// Set input and output for the registration filter
		registration->SetInput(target, source);
		registration->SetOutput(mffd_ed_b[t]);

		// Make an initial Guess for the parameters.
		registration->GuessParameter();
		// Overrride with any the user has set.
		if (parin_name != NULL) {
			registration->irtkImageRegistration2::Read(parin_name);
		}

		// Run registration filter
		registration->Run();

		if (registration->GetMFFDMode()) {
			mffd_ed_b[t]->CombineLocalTransformation();
		}

		if (debug) {
			char buffer[255];
			sprintf(buffer, "debug_ed_b_%.2d.dof.gz", t);
			mffd_ed_b[t]->irtkTransformation::Write(buffer);
		}

		delete registration;
		delete target;
		delete source;
	}
	/// Weighted Combine forward backward
	for(i = 0; i < cine.GetT(); i++){

		for(j = 0; j < mffd_ed_f[i]->GetLocalTransformation(0)->NumberOfDOFs(); j++){
			double newdof, fdof, bdof;
			fdof = mffd_ed_f[i]->GetLocalTransformation(0)->Get(j);
			bdof = mffd_ed_b[i]->GetLocalTransformation(0)->Get(j);
			newdof = ((cine.GetT() - i - 1)*fdof + i*bdof)/(cine.GetT() - 1);
			mffd_ed_f[i]->GetLocalTransformation(0)->Put(j,newdof);
		}

		if (debug) {
			char buffer[255];
			sprintf(buffer, "debug_ed_%.2d.dof.gz", i);
			mffd_ed_f[i]->irtkTransformation::Write(buffer);
		}
	}
	/// Third, propergate from ES to all
	/// Forward
	// create transformations and registrations
	irtkMultiLevelFreeFormTransformation **mffd_es_f = new irtkMultiLevelFreeFormTransformation*[cine.GetT()];
	for(i = 0; i < cine.GetT(); i++){
		oldt = t;
		t = esphase+i;
		if(t > cine.GetT() - 1){
			t -= cine.GetT();
		}

		if(i > 0){
			mffd_es_f[t] = new irtkMultiLevelFreeFormTransformation(*mffd_es_f[oldt]);
		}else{
			mffd_es_f[t] = new irtkMultiLevelFreeFormTransformation;
		}

		irtkImageFreeFormRegistration2 *registration = new irtkImageFreeFormRegistration2;

		// Combine images
		irtkImageAttributes attr = cine.GetImageAttributes();
		attr._t = 1;
		irtkGreyImage *target = new irtkGreyImage(attr);
		irtkGreyImage *source = new irtkGreyImage(attr);
		for (z = 0; z < target->GetZ(); z++) {
			for (y = 0; y < target->GetY(); y++) {
				for (x = 0; x < target->GetX(); x++) {
					target->Put(x, y, z, 0, cine.Get(x, y, z, esphase));
					source->Put(x, y, z, 0, cine.Get(x, y, z, t));
				}
			}
		}

		// Set input and output for the registration filter
		registration->SetInput(target, source);
		registration->SetOutput(mffd_es_f[t]);

		// Make an initial Guess for the parameters.
		registration->GuessParameter();
		// Overrride with any the user has set.
		if (parin_name != NULL) {
			registration->irtkImageRegistration2::Read(parin_name);
		}

		// Run registration filter
		registration->Run();

		if (registration->GetMFFDMode()) {
			mffd_es_f[t]->CombineLocalTransformation();
		}

		if (debug) {
			char buffer[255];
			sprintf(buffer, "debug_es_f_%.2d.dof.gz", t);
			mffd_es_f[t]->irtkTransformation::Write(buffer);
		}

		delete registration;
		delete target;
		delete source;
	}
	/// Backward
	irtkMultiLevelFreeFormTransformation **mffd_es_b = new irtkMultiLevelFreeFormTransformation*[cine.GetT()];
	for(i = 0; i < cine.GetT(); i++){
		oldt = t;
		t = esphase-i;
		if(t < 0){
			t += cine.GetT();
		}

		if(i > 0){
			mffd_es_b[t] = new irtkMultiLevelFreeFormTransformation(*mffd_es_b[oldt]);
		}else{
			mffd_es_b[t] = new irtkMultiLevelFreeFormTransformation;
		}

		irtkImageFreeFormRegistration2 *registration = new irtkImageFreeFormRegistration2;

		// Combine images
		irtkImageAttributes attr = cine.GetImageAttributes();
		attr._t = 1;
		irtkGreyImage *target = new irtkGreyImage(attr);
		irtkGreyImage *source = new irtkGreyImage(attr);
		for (z = 0; z < target->GetZ(); z++) {
			for (y = 0; y < target->GetY(); y++) {
				for (x = 0; x < target->GetX(); x++) {
					target->Put(x, y, z, 0, cine.Get(x, y, z, esphase));
					source->Put(x, y, z, 0, cine.Get(x, y, z, t));
				}
			}
		}

		// Set input and output for the registration filter
		registration->SetInput(target, source);
		registration->SetOutput(mffd_es_b[t]);

		// Make an initial Guess for the parameters.
		registration->GuessParameter();
		// Overrride with any the user has set.
		if (parin_name != NULL) {
			registration->irtkImageRegistration2::Read(parin_name);
		}

		// Run registration filter
		registration->Run();

		if (registration->GetMFFDMode()) {
			mffd_es_b[t]->CombineLocalTransformation();
		}

		if (debug) {
			char buffer[255];
			sprintf(buffer, "debug_es_b_%.2d.dof.gz", t);
			mffd_es_b[t]->irtkTransformation::Write(buffer);
		}

		delete registration;
		delete target;
		delete source;
	}
	/// Weighted Combined forward backward
	for(i = 0; i < cine.GetT(); i++){

		t = i - esphase;
		if(t < 0){
			t += cine.GetT();
		}

		for(j = 0; j < mffd_es_f[i]->GetLocalTransformation(0)->NumberOfDOFs(); j++){
			double newdof, fdof, bdof;
			fdof = mffd_es_f[i]->GetLocalTransformation(0)->Get(j);
			bdof = mffd_es_b[i]->GetLocalTransformation(0)->Get(j);

			newdof = ((cine.GetT() - t - 1)*fdof + t*bdof)/(cine.GetT() - 1);
			mffd_es_f[i]->GetLocalTransformation(0)->Put(j,newdof);
		}

		if (debug) {
			char buffer[255];
			sprintf(buffer, "debug_es_%.2d.dof.gz", i);
			mffd_es_f[i]->irtkTransformation::Write(buffer);
		}
	}
	/// Fourth patch based segmentation from propergated ED and ES
	for(i = 0; i < cine.GetT(); i++){
		irtkImageFunction *interpolator = 
			new irtkNearestNeighborInterpolateImageFunction;
		irtkImageFunction *interpolator2 = 
			new irtkLinearInterpolateImageFunction;
		// transform segmentation from ED to current
		irtkImageTransformation *imagetransformation =
			new irtkImageTransformation;

		imagetransformation->SetInput (&segmentation_ed, mffd_ed_f[i]);
		imagetransformation->SetOutput(&propergated_segmentation_ed);
		imagetransformation->PutInterpolator(interpolator);
		imagetransformation->InvertOn();

		// Transform image
		imagetransformation->Run();

		irtkRealImage tmp(blured.GetFrame(0));

		imagetransformation->SetInput (&tmp, mffd_ed_f[i]);
		imagetransformation->SetOutput(&propergated_image_ed);
		imagetransformation->PutInterpolator(interpolator2);
		imagetransformation->InvertOn();

		// Transform image
		imagetransformation->Run();

		// transform segmentation from ES to current

		imagetransformation->SetInput (&segmentation_es, mffd_es_f[i]);
		imagetransformation->SetOutput(&propergated_segmentation_es);
		imagetransformation->PutInterpolator(interpolator);
		imagetransformation->InvertOn();

		// Transform image
		imagetransformation->Run();

		irtkRealImage tmp2(blured.GetFrame(esphase));

		imagetransformation->SetInput (&tmp2, mffd_es_f[i]);
		imagetransformation->SetOutput(&propergated_image_es);
		imagetransformation->PutInterpolator(interpolator2);
		imagetransformation->InvertOn();

		// Transform image
		imagetransformation->Run();

		tmp_image = blured.GetFrame(i);
		// run patch based segmentation
		irtkRealImage **atlases = new irtkRealImage*[2];
		atlases[0] = &propergated_image_ed;
		atlases[1] = &propergated_image_es;
		irtkGreyImage **labels = new irtkGreyImage*[2];
		labels[0] = &propergated_segmentation_ed;
		labels[1] = &propergated_segmentation_es;

		if(debug){
			tmp_image.Write("patch_image.nii.gz");
			atlases[0]->Write("patch_atlas0.nii.gz");
			atlases[1]->Write("patch_atlas1.nii.gz");
			labels[0]->Write("patch_label0.nii.gz");
			labels[1]->Write("patch_label1.nii.gz");
		}

		irtkPatchBasedSegmentation patchBased(tmp_image, atlases, labels, 2, 2, 2);
		patchBased.SetPadding(-1);
		patchBased.Run();
		patchBased.GetCardiacSegmentation(cine.GetT()/2 - min(abs(i),abs(cine.GetT() - i)), 
			cine.GetT()/2 - min(abs(i - esphase),abs(cine.GetT() - i + esphase)));
		tmp_segmentation = patchBased.GetSegmentation();

		if(debug){
			tmp_segmentation.Write("patch_output.nii.gz");
		}

		for (z = 0; z < cine.GetZ(); z++) {
			for (y = 0; y < cine.GetY(); y++) {
				for (x = 0; x < cine.GetX(); x++) {
					output.Put(x, y, z, i, tmp_segmentation.Get(x, y, z));
				}
			}
		}

		delete imagetransformation;
		delete interpolator;
		delete interpolator2;
		delete []labels;
		delete []atlases;
	}

	output.Write(output_name);
}
