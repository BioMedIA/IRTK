/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

// irtk includes
#include <irtkImage.h>

#include <irtkNoise.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: noise [in] [out] [noise] <parameters>\n";
  cerr << "where [noise] is one of the following:\n";
  cerr << "\t -uniform             Uniform noise\n";
  cerr << "\t -gaussian            Gaussian white noise\n";
  cerr << "\t -rician              Rician noise\n";
  cerr << "and where <parameters> is one of the following:\n";
  cerr << "\t -Tp                  Padding value (for all noise types)\n";
  cerr << "\t -amplitude           Amplitude/STD of uniform or Rician noise\n";
  cerr << "\t -percent             Percentage for highest signal as Amplitude/STD of uniform or Rician noise\n";
  cerr << "\t -mean                Mean  of Gaussian noise\n";
  cerr << "\t -sigma               Sigma of Gaussian noise\n";
  cerr << "\t -minval              Min value of Gaussian noise\n";
  cerr << "\t -maxval              Max value of Gaussian noise\n\n";
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  irtkNoise<irtkGreyPixel>                    *noise          = NULL;
  irtkUniformNoiseWithPadding<irtkGreyPixel>  *uniform_noise  = NULL;
  irtkGaussianNoiseWithPadding<irtkGreyPixel> *gaussian_noise = NULL;
  irtkRicianNoiseWithPadding<irtkGreyPixel>   *rician_noise   = NULL;
  irtkGreyImage input;

  if (argc < 4) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  input.Read(input_name);

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-uniform") == 0)) {
      argc--;
      argv++;
      if (noise == NULL) {
        noise = uniform_noise = new irtkUniformNoiseWithPadding<irtkGreyPixel>;
      } else {
        cerr << "More than one type of noise specified" << endl;
        usage();
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-gaussian") == 0)) {
      argc--;
      argv++;
      if (noise == NULL) {
        noise = gaussian_noise = new irtkGaussianNoiseWithPadding<irtkGreyPixel>;
      } else {
        cerr << "More than one type of noise specified" << endl;
        usage();
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-rician") == 0)) {
      argc--;
      argv++;
      if (noise == NULL) {
#ifdef OLD_RICIAN
        cerr << "Using old version of Rician noise\n";
#else
        cerr << "Using corrected version of Rician noise\n";
#endif
        noise = rician_noise = new irtkRicianNoiseWithPadding<irtkGreyPixel>;
      } else {
        cerr << "More than one type of noise specified" << endl;
        usage();
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-amplitude") == 0)) {
      argc--;
      argv++;
      if (noise != NULL) {
        noise->SetAmplitude(atof(argv[1]));
      } else {
        cerr << "No noise type specified" << endl;
        usage();
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-percent") == 0)) {
      argc--;
      argv++;
      irtkGreyPixel min, max;
      input.GetMinMax(&min, &max);
      if (noise != NULL) {
        cerr << "Setting noise amplitude to " << atof(argv[1])*double(max) << endl;
        noise->SetAmplitude(atof(argv[1])*double(max));
      } else {
        cerr << "No noise type specified" << endl;
        usage();
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mean") == 0)) {
      argc--;
      argv++;
      if (gaussian_noise != NULL) {
        gaussian_noise->SetMean(atof(argv[1]));
      } else {
        cerr << "No gaussian noise type specified for option -mean" << endl;
        usage();
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)) {
      argc--;
      argv++;
      if (gaussian_noise != NULL) {
        gaussian_noise->SetSigma(atof(argv[1]));
      } else {
        cerr << "No gaussian noise type specified for option -sigma" << endl;
        usage();
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-minval") == 0)) {
      argc--;
      argv++;
      if (gaussian_noise != NULL) {
        gaussian_noise->SetMinVal(atoi(argv[1]));
      } else {
        cerr << "No gaussian noise type specified for option -minval" << endl;
        usage();
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-maxval") == 0)) {
      argc--;
      argv++;
      if (gaussian_noise != NULL) {
        gaussian_noise->SetMaxVal(atoi(argv[1]));
      } else {
        cerr << "No gaussian noise type specified for option -maxval" << endl;
        usage();
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      if (noise != NULL) {
        if (uniform_noise != NULL) {
          uniform_noise->SetPaddingValue(atoi(argv[1]));
        } else if (gaussian_noise != NULL) {
          gaussian_noise->SetPaddingValue(atoi(argv[1]));
        } else if (rician_noise != NULL) {
          rician_noise->SetPaddingValue(atoi(argv[1]));
        }
      } else {
        cerr << "No noise type specified" << endl;
        usage();
      }
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Couldn't parse args" << endl;
      usage();
    }
  }

  // Noise
  noise->SetInput(&input);
  noise->SetOutput(&input);
  noise->Run();

  // Write output
  input.Write(output_name);

  return 0;
}

