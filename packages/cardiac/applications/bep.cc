#include <irtkBep.h>

#ifdef HAS_VTK

char *bep_name = NULL, *data_name = NULL, *output_name;

void usage()
{
	cerr << "Usage: bep [Number of frames] [Number of segments] [Bep Surface] [Data Surface] [Output]" << endl;
	cerr << "multiple frame, data surface is prefix name, actual data name prefix+%2d" << endl;
	exit(1);
}

int main( int argc, char** argv )
{
	int i, numberofframes, numberofsegments;
	irtkBep bep;

	if (argc < 6) {
		usage();
	}

	// Parse filenames
	numberofframes = atoi(argv[1]);
	argv++;
	argc--;
    numberofsegments = atoi(argv[1]);
    argv++;
    argc--;
	bep_name = argv[1];
	argv++;
	argc--;
	data_name = argv[1];
	argv++;
	argc--;
	output_name = argv[1];
	argv++;
	argc--;

    // bepsurface pipeline
    cout << "Reading bep surface ... " << bep_name << endl;
    vtkPolyDataReader *bepsurface_reader = vtkPolyDataReader::New();
    bepsurface_reader->SetFileName(bep_name);
    bepsurface_reader->Modified();
    bepsurface_reader->Update();
    vtkPolyData *bepsurface = vtkPolyData::New();
    bepsurface = bepsurface_reader->GetOutput();
    bepsurface->Update();

	//go
	if(numberofframes == 1){

		// datasurface pipeline
		cout << "Reading data surface ... " << data_name << endl;
		vtkPolyDataReader *datasurface_reader = vtkPolyDataReader::New();
		datasurface_reader->SetFileName(data_name);
		datasurface_reader->Modified();
		datasurface_reader->Update();
		vtkPolyData *datasurface = vtkPolyData::New();
		datasurface = datasurface_reader->GetOutput();
		datasurface->Update();

		bep.SetInput(bepsurface,datasurface,numberofsegments);
		bep.SetOutput(output_name);
		bep.Initialize();
		bep.Bullseyeplot();
		bep.Finalize();

		datasurface_reader->Delete();
	}else{
		for(i = 0; i < numberofframes; i++){

            char buffer[255];
            sprintf(buffer, "%s%.2d.vtk", data_name, i);

            // datasurface pipeline
            cout << "Reading data surface ... " << buffer << endl;
            vtkPolyDataReader *datasurface_reader = vtkPolyDataReader::New();
            datasurface_reader->SetFileName(buffer);
            datasurface_reader->Modified();
            datasurface_reader->Update();
            vtkPolyData *datasurface = vtkPolyData::New();
            datasurface->DeepCopy(datasurface_reader->GetOutput());
            datasurface->Update();

            bep.SetInput(bepsurface,datasurface,numberofsegments);
            bep.SetOutput(output_name);
            bep.Initialize();
            bep.Bullseyeplot();
            bep.Finalize();

            datasurface->Delete();
            datasurface_reader->Delete();
		}
	}
    bepsurface_reader->Delete();
}
#else
#include <irtkCommon.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif