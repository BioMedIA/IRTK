


//#include <string>
#include <irtkRegistration.h>



using namespace std;

void usage ()
{
	cerr << "Usage: ComputeTFFDError [solutionFile] [Nt] [errorFiles_n]" << endl;
	exit(1);
}


int main(int argc, char **argv)
{
	if (argc < 4) {
		usage();
	  }

	char * solution_name = NULL;
	char ** errorFile_name = NULL;
	int Nt, Np, NpT, i;

	double error;
    double errorMean = 0;
    double errorSTD = 0;
    double errorMin = 100;
    double errorMax = 0;
    double * errorMeanT = NULL;
    double * errorSTDT = NULL;
    double * errorMinT = NULL;
    double * errorMaxT = NULL;

	cout<<"reading input ... "; cout.flush();
	solution_name = argv[1];
	argc--;
	argv++;
	Nt = atoi(argv[1]);
	argc--;
	argv++;
	errorFile_name = new char *[Nt];
	for (i = 0; i < Nt; i++) {
		errorFile_name[i] = argv[1];
		argc--;
		argv++;
	}
	errorMeanT = new double [Nt];
	errorSTDT = new double [Nt];
	errorMinT = new double [Nt];
	errorMaxT = new double [Nt];
	cout<<"done"<<endl;

	cout<<"compute error metrics ... "; cout.flush();
	Np = 0;
	// compute mean, min, max:
	for (i = 0; i < Nt; i++) {
	  NpT = 0;
	  errorMinT[i] = 100;
	  errorMaxT[i] = 0;
	  ifstream textFileReader;
	  textFileReader.open(errorFile_name[i]);
	  if (textFileReader.is_open()) {
		while (!textFileReader.eof()) {
		  textFileReader >> ws;
		  textFileReader >> error;
		  textFileReader >> ws;

		  errorMean += error;
		  errorMin = (errorMin > error) ? error : errorMin;
		  errorMax = (errorMax < error) ? error : errorMax;
		  errorMeanT[i] += error;
		  errorMinT[i] = (errorMinT[i] > error) ? error : errorMinT[i];
		  errorMaxT[i] = (errorMaxT[i] < error) ? error : errorMaxT[i];
		  Np++;
		  NpT++;
		}
	  } else {
		cerr<<"coudn't open time file "<<errorFile_name[i]<<endl;
		exit(1);
	  }
	  errorMeanT[i] /= NpT;
	  textFileReader.close();
	}
	errorMean /= Np;
	// compute std:
	for (i = 0; i < Nt; i++) {
	  NpT = 0;
	  ifstream textFileReader;
	  textFileReader.open(errorFile_name[i]);
	  if (textFileReader.is_open()) {
		while (!textFileReader.eof()) {
		  textFileReader >> ws;
		  textFileReader >> error;
		  textFileReader >> ws;

		  errorSTD += (error - errorMean) * (error - errorMean);
		  errorSTDT[i] += (error - errorMeanT[i]) * (error - errorMeanT[i]);
		  NpT++;
		}
	  } else {
		cerr<<"coudn't open time file "<<errorFile_name[i]<<endl;
		exit(1);
	  }
	  errorSTDT[i] /= NpT;
	  errorSTDT[i] = sqrt(errorSTDT[i]);
	  textFileReader.close();
	}
	errorSTD /= Np;
	errorSTD = sqrt(errorSTD);
	cout<<"done"<<endl;

    cout << "error all mean: " << errorMean << endl;
    cout << "error all std:  " << errorSTD  << endl;
    cout << "error all min:  " << errorMin  << endl;
    cout << "error all max:  " << errorMax  << endl;

	fstream to(solution_name, ios::out | ios::app);
	if (!to) {
		cerr << "Writing solution file: Can't open file " << solution_name;
		exit(1);
	}
	for (i = 0; i < Nt; i++) {
	  to << "error for frame" << i << "\n";
	  to << "error mean: " << errorMeanT[i] << "\n";
	  to << "error std:  " << errorSTDT[i]  << "\n";
	  to << "error min:  " << errorMinT[i]  << "\n";
	  to << "error max:  " << errorMaxT[i]  << "\n";
	  to << "\n";
	}
	to << "---------------------------------------\n";
	to << "error all mean: " << errorMean << "\n";
	to << "error all std:  " << errorSTD  << "\n";
	to << "error all min:  " << errorMin  << "\n";
	to << "error all max:  " << errorMax  << "\n";
	to.close();

	return 0;
}









