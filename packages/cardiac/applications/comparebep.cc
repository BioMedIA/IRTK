#include <irtkCommon.h>

void usage()
{
  cerr << "Usage: comparebep [ResultInput1] [ResultInput2] [ResultOutput] [MaxOutput] [numberofframes]" << endl;
  exit(1);
}

int main( int argc, char** argv )
{
	char *inputfilename1 = NULL; char *inputfilename2 = NULL;
	char *resultfilename = NULL; char *maxfilename = NULL;
	int number = 20,t,i;

	double input1,input2,result,max = 0;

	if (argc<5) usage();

	inputfilename1 = argv[1];
	argc--;
	argv++;
	inputfilename2 = argv[1];
	argc--;
	argv++;
	resultfilename = argv[1];
	argc--;
	argv++;
	maxfilename = argv[1];
	argc--;
	argv++;
	if(argc > 1)
	number = atoi(argv[1]);
	argc--;
	argv++;

	//openfiles
	remove(maxfilename);
	remove(resultfilename);

	ifstream inputin1(inputfilename1, ios::in);
	ifstream inputin2(inputfilename2, ios::in);
	ofstream maxout(maxfilename,ios::app);
	ofstream resultout(resultfilename,ios::app);
	

	for(t=0;t<number-1;t++){
		  for(i=0; i<17; i++){
			  inputin1 >> input1;
			  inputin2 >> input2;
			  result = input2 - input1;
			  resultout<<result<<" ";
			  if(abs(result) > max) max = abs(result);
		  }
		  resultout<<endl;
	}
	maxout<<max<<" ";

	for(i=0; i<17; i++){
		inputin1 >> input1;
		inputin2 >> input2;
		result = input2 - input1;
		resultout<<result<<" ";
		if(abs(result) > max) max = abs(result);
	}
	resultout<<endl;
	resultout.close();
	maxout<<max<<endl;
	maxout.close();
}