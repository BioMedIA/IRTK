
#include <irtkRegistration.h>
#include <irtkSpectralClustering.h>

irtkMatrix read_csv(string csvFilename1, string csvFilename2);
void writeAtlases(string filename, int * atlases, int nrAtlases);
void writeAtlasNames(string filename, int * atlases, int nrAtlases, vector<string> atlasNames);
void usage()
{
	cerr << "usage LEAPatlas_selection [Atlas Similarities] [Target Similarities] [N selected atlases] [output]" << endl;
	cerr << "-atlasNames [list of atlas names]. With this option, atlas names are written instead of numbers." << endl;

		
}
int main(int argc, char **argv)
{
	if(argc < 4){
		usage();
		exit(1);
	}

	string atlassims = argv[1];
	argc--;
	argv++;
	string targetsims = argv[1];
	argc--;
	argv++;
	int nrAtlases = atoi(argv[1]);
	argc--;
	argv++;
	string outname = argv[1];
	argc--;
	argv++;


	bool ok;
	bool useAtlasNames = false;
	vector<string> atlasNames = vector<string>();
	while (argc > 1){
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-atlasNames") == 0)){
		  argc--;
		  argv++;

		  string line;
		  string imagefilename = argv[1];
		  argc--;
		  argv++;
		  cout << "Using atlas names from " << imagefilename << endl;
		  ifstream from(imagefilename.c_str());
		  if (!from) {
			cerr << "Can't open file "
			<< endl;
			exit(1);
		  }
		  if (from.is_open())	{
			while (! from.eof() ) {
			  getline (from,line);
			  string fileName = line;
			  atlasNames.push_back(fileName);
			}
			from.close();
		  }
		  useAtlasNames = true;
		  ok = true;
		}
		if (ok == false){
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}
	cout << "Read similarity matrices" << endl;
	irtkMatrix w;
	try{
		w = read_csv(atlassims, targetsims);
	}
	catch (int e){
		if(e == 1){
			cerr << "File not found exception" << endl;
		}
		return -1;
	}
	irtkSpectralClustering le(w.Rows(), 2);
	cout << "Do manifold embedding" << endl;
	le.Initialize(w);
	le.DoSpectralEmbedding();
	int selected_atlases[nrAtlases];
	cout << "Select " << nrAtlases << " closest atlases (IDs in [0...nAtlases-1]): " << endl;
	le.GetNeighbours(w.Rows()-1, nrAtlases, selected_atlases);

	if(useAtlasNames)
		writeAtlasNames(outname, selected_atlases, nrAtlases, atlasNames);
	else
		writeAtlases(outname, selected_atlases, nrAtlases);
}

irtkMatrix read_csv(string csvFilename1, string csvFilename2){

	ifstream input;
	input.open(csvFilename1.c_str());
	string value;
	int rows = 0;
	if(! input){
		cerr << "File " << csvFilename1 << " does not exist" << endl;
		throw 1;
	}
	while(input.good()){
		getline(input, value, '\n');
		rows++;

	}
	rows--;
//	if(rows != _nrSubjects){
//		cerr << "Dimensions in .csv file do not agree with nrSubjects" << endl;
//	}
	irtkMatrix _w;
	_w.Initialize(rows+1, rows+1);

	ifstream input2;
	input2.open(csvFilename1.c_str());
	int rowctr = 0;
	int colctr = 0;
	while(input2.good()){
		getline(input2, value, '\n');
		int pos = value.find(",");
		int posOld = -1;
		while(colctr < rows && rowctr < rows){
			string stringdata = value.substr(posOld+1,pos-posOld-1);
			posOld = pos;
			pos = value.find(",", posOld+1);
			double data = atof(stringdata.c_str());
			_w(rowctr, colctr) = data;
			colctr++;
		}
		colctr = 0;
		rowctr ++;
//		cout << rowctr << endl;
	}
//	cout << "a" << endl;
	ifstream input3;
	input3.open(csvFilename2.c_str());
	if(! input3){
			cerr << "File " << csvFilename2 << " does not exist" << endl;
			throw 1;
		}
	rowctr = 0;
	while(input3.good() && rowctr < rows){
		getline(input3, value, '\n');
		double data = atof(value.c_str());
		_w(rowctr, rows) = data;
		_w(rows, rowctr) = data;
		rowctr++;
//		cout << data << endl;
	}
	if(rowctr != rows){
		cerr << "No of target similarities (" << rowctr << ") does not agree with No of atlas similarities(" << rows << ")"<< endl;

	}
	_w(rows, rows) = 2;


	return _w;

}

void writeAtlasNames(string filename, int * atlases, int nrAtlases, vector<string> atlasNames){
	ofstream output;
	cout << "Write atlas names to '" << filename << "'" << endl;
	output.open(filename.c_str());
	for(int r = 0; r < nrAtlases; r++){
		output << atlasNames[*atlases];
		output << endl;
		atlases++;

	}
	output.close();
}

void writeAtlases(string filename, int * atlases, int nrAtlases){
	ofstream output;
	cout << "Write atlas IDs to '" << filename << "'" << endl;
	output.open(filename.c_str());
	for(int r = 0; r < nrAtlases; r++){


		output << *atlases;
	//	cout << *atlases << endl;
		if(r != nrAtlases-1)
			output << ","; //endl;
		atlases++;

	}

	output << endl;

	output.close();
}
