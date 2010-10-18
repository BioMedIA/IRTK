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

#include <irtkHistogram.h>

#include <irtkTransformation.h>

#include <nr.h>

// Definition of available adaptivity measures
typedef enum { Entropy, Variance } AdaptivityMeasure;

// Definition of available adaptivity selection methods
typedef enum { PercentOfDOFs,
               PercentOfAdaptivity
             } AdaptivitySelection;

// Input transformation
char *dofin_name  = NULL;

// Output transformation
char *dofout_name = NULL;

// Image for adavptivity calculation
char *target_name = NULL;

irtkGreyImage *target;
irtkFreeFormTransformation3D *ffd;
irtkMultiLevelFreeFormTransformation *mffd;

void usage()
{
  cerr << "Usage: ffdedit [target] [dofin] [dofout] [options]" << endl;
  cerr << "where [options] are one or more of the following: " << endl;
  cerr << "-mode     <n>" << endl;
  cerr << "-level    <n>" << endl;
  cerr << "-percent  <n>" << endl;
  cerr << "-measure  <n>" << endl;
  exit(1);
}

void FFDEdit1(AdaptivityMeasure measure, double value)
{
  float *arr, *brr;
  double adaptivity;
  irtkGreyPixel min, max, *ptr2target;
  int index, x, y, z, x1, y1, z1, x2, y2, z2, n;

  // Figure out number of control points
  n = ffd->NumberOfDOFs()/3;

  // Allocate
  arr = new float [n];
  brr = new float [n];

  // Set up 1D histogram
  target->GetMinMax(&min, &max);
  irtkHistogram histo(max-min+1);
  histo.PutMin(min - 0.5);
  histo.PutMax(max + 0.5);

  // Get lattice dimensions and loop over lattice
  for (x = 0; x < ffd->GetX(); x++) {
    for (y = 0; y < ffd->GetY(); y++) {
      for (z = 0; z < ffd->GetZ(); z++) {

        // Get index and point from lattice coordinates
        index = ffd->LatticeToIndex(x, y, z);
        ffd->BoundingBox(target, index, x1, y1, z1, x2, y2, z2);

        // Reset histogram and fill
        histo.Reset();
        for (int k = z1; k <= z2; k++) {
          for (int j = y1; j <= y2; j++) {

            // Get pointer to target image
            ptr2target = target->GetPointerToVoxels(x1, j, k);

            for (int i = x1; i <= x2; i++) {
              histo.AddSample(*ptr2target);
              ptr2target++;
            }
          }
        }

        // Compute adaptivity measure
        switch (measure) {
        case Variance:
          adaptivity = histo.Variance();
          break;
        case Entropy:
          adaptivity = histo.Entropy();
          break;
        default:
          cerr << "Unknown adaptivity measure" << endl;
          exit(1);
        }

        // Store in array
        arr[index] = adaptivity;
        brr[index] = index;
      }
    }
  }

  // Quicksort
  sort2(n, arr-1, brr-1);

  // Set top % to active, bottom % to passive
  index = round(value*(double)n);
  for (int i=(int)n; i>index; i--) {
    ffd->IndexToLattice(int(brr[i-1]), x, y, z);
    ffd->PutStatus(x, y, z, _Active);
  }
  for (int i=index; i>0; i--) {
    ffd->IndexToLattice(int(brr[i-1]), x, y, z);
    ffd->PutStatus(x, y, z, _Passive);
  }

  // Be good
  delete [] arr;
  delete [] brr;
}

void FFDEdit2(AdaptivityMeasure measure, double value)
{
  double adaptivity;
  irtkGreyPixel min, max, *ptr2target;
  int index, x, y, z, x1, y1, z1, x2, y2, z2;

  // Set up 1D histogram
  target->GetMinMax(&min, &max);
  irtkHistogram histo(max-min+1);
  histo.PutMin(min - 0.5);
  histo.PutMax(max + 0.5);

  // Get lattice dimensions and loop over lattice
  ptr2target = target->GetPointerToVoxels();
  for (int z=0; z<target->GetZ(); z++) {
    for (int y=0; y<target->GetY(); y++) {
      for (int x=0; x<target->GetX(); x++) {
        histo.AddSample(*ptr2target);
        ptr2target++;
      }
    }
  }

  // Compute adaptivity measure
  switch (measure) {
  case Variance:
    value = value * histo.Variance();
    break;
  case Entropy:
    value = value * histo.Entropy();
    break;
  default:
    cerr << "Unknown adaptivity measure" << endl;
    exit(1);
  }

  for (x = 0; x < ffd->GetX(); x++) {
    for (y = 0; y < ffd->GetY(); y++) {
      for (z = 0; z < ffd->GetZ(); z++) {

        // Get index and point from lattice coordinates
        index = ffd->LatticeToIndex(x, y, z);
        ffd->BoundingBox(target, index, x1, y1, z1, x2, y2, z2);

        // Reset histogram and fill
        histo.Reset();
        for (int k = z1; k <= z2; k++) {
          for (int j = y1; j <= y2; j++) {

            // Get pointer to target image
            ptr2target = target->GetPointerToVoxels(x1, j, k);

            for (int i = x1; i <= x2; i++) {
              histo.AddSample(*ptr2target);
              ptr2target++;
            }
          }
        }

        // Compute adaptivity measure
        switch (measure) {
        case Variance:
          adaptivity = histo.Variance();
          break;
        case Entropy:
          adaptivity = histo.Entropy();
          break;
        default:
          cerr << "Unknown adaptivity measure" << endl;
          exit(1);
        }

        // Set status according to adaptivity in region
        if (adaptivity < value) {
          ffd->PutStatus(x, y, z, _Passive);
        } else {
          ffd->PutStatus(x, y, z, _Active);
        }
      }
    }
  }
}

int main(int argc, char **argv)
{
  double percentage;
  int i, x, y, z, ok, active, passive, level;
  AdaptivityMeasure measure;
  AdaptivitySelection selection;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse file names
  target_name = argv[1];
  argc--;
  argv++;
  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read target
  cout << "Reading image ... "; cout.flush();
  target = new irtkGreyImage;
  target->Read(target_name);
  cout << "done" << endl;

  // Read transformation
  mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);

  // Default level is the last level
  level = mffd->NumberOfLevels()-1;

  // Default measure is entropy
  measure = Entropy;

  // Default selection is percent of DOFs
  selection = PercentOfDOFs;

  // Default percenttage is 50%
  percentage = 0.5;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-measure") == 0)) {
      argc--;
      argv++;
      if (strcmp(argv[1], "entropy") == 0) {
        measure = Entropy;
      } else {
        if (strcmp(argv[1], "variance") == 0) {
          measure = Variance;
        } else {
          cerr << "Unknown measure: " << argv[1] << endl;
          exit(1);
        }
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mode") == 0)) {
      argc--;
      argv++;
      if (strcmp(argv[1], "PercentOfDOFs") == 0) {
        selection = PercentOfDOFs;
      } else {
        if (strcmp(argv[1], "PercentOfAdaptivity") == 0) {
          selection = PercentOfAdaptivity;
        } else {
          cerr << "Unknown mode: " << argv[1] << endl;
          exit(1);
        }
      }
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-percent") == 0)) {
      argc--;
      argv++;
      percentage = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-level") == 0)) {
      argc--;
      argv++;
      level = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can't parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Extract current transformation level
  ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(level));

  if (ffd == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

  // Edit current transformation level
  cout << "Editing transformation level " << level << " ... ";
  cout.flush();
  switch (selection) {
  case PercentOfDOFs:
    FFDEdit1(measure, percentage);
    break;
  case PercentOfAdaptivity:
    FFDEdit2(measure, percentage);
    break;
  }
  cout << "done" << endl;

  // Calculate number of active and passive control points
  i = 0;
  active  = 0;
  passive = 0;
  for (x = 0; x < ffd->GetX(); x++) {
    for (y = 0; y < ffd->GetY(); y++) {
      for (z = 0; z < ffd->GetZ(); z++) {
        _Status sx, sy, sz;
        ffd->GetStatus(x, y, z, sx, sy, sz);
        if (sx == _Active) {
          active++;
        } else {
          passive++;
        }
        if (sy == _Active) {
          active++;
        } else {
          passive++;
        }
        if (sz == _Active) {
          active++;
        } else {
          passive++;
        }
        i++;
      }
    }
  }
  cout << "Transformation level " << level << " has " << endl;
  cout << active << " active control points  \t ("
       << active*100.0/(active+passive) << "%)" << endl;
  cout << passive << " passive control points \t ("
       << passive*100.0/(active+passive) << "%)" << endl;

  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}


