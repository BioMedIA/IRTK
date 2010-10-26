/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  pa100 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  $
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

char *imageA_name = NULL;
char *imageB_name = NULL;
char *labelFile   = NULL;

#define MAX_LABELS 1000

void usage()
{
  cerr << " labelStats [imageA] [imageB] <-options>" << endl;
  cerr << " " << endl;
  cerr << " Print for each label data on number and overlap in the images."  << endl;
  cerr << " Images contain structural labels (max 1000 labels)." << endl;
  cerr << " Images must have same dimensions."  << endl;
  cerr << " Output format is comma separated and one row per label:" << endl;
  cerr << " " << endl;
  cerr << "          Label Number , n(A) , n(B) , n(A^B) , OR , Dice" << endl;
  cerr << " " << endl;
  cerr << " Where OR is Tanimoto overlap and Dice is the Dice coefficient or overlap" << endl;
  cerr << " " << endl;
  cerr << " options:" << endl;
  cerr << " " << endl;
  cerr << " -summary : print only the mean overlaps over all labels." << endl;
  cerr << " -file    : File containg a restricted list of labels (separated by whitespace)." << endl;
  cerr << " -false   : print false positives / negatives, format becomes:" << endl;
  cerr << "            Label,n(A),n(B),n(A^B),n(A\\B),n(B\\A),OR,Dice" << endl;
  cerr << " -q       : Quiet summary, csv output : 'OR,Dice' " << endl;
  cerr << " -diceRow : Comma separated Dice values for all labels in single a row (no labels printed)." << endl;
  cerr << " -orRow   : Comma separated Tanimoto overlaps for all labels in single a row (no labels printed)." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, nVoxels;
  irtkGreyPixel *pA, *pB, max, min, tmp;
  int jointLabels[MAX_LABELS][MAX_LABELS];
  int marginalA[MAX_LABELS];
  int marginalB[MAX_LABELS];
  int includeLabel[MAX_LABELS];

  bool printSummary = false;
  bool quiet = false;
  bool printFalseValues = false;
  bool ok;
  int temp;
  bool diceRow = false;
  bool orRow = false;

  double overlap, dice;
  double sumOverlap, sumDice;
  int count;

  // Read arguments
  if (argc < 3){
    usage();
  }

  imageA_name = argv[1];
  argc--;
  argv++;
  imageB_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-summary") == 0)){
      argc--;
      argv++;
      printSummary = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-diceRow") == 0)){
      argc--;
      argv++;
      diceRow = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-orRow") == 0)){
      argc--;
      argv++;
      orRow = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-file") == 0)){
      argc--;
      argv++;
      labelFile = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-false") == 0)){
      argc--;
      argv++;
      printFalseValues = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  if (quiet == true){
    // Force summary.
    printSummary = true;
  }


  irtkGreyImage *imageA, *imageB;
  imageA = new irtkGreyImage(imageA_name);
  imageB = new irtkGreyImage(imageB_name);

  if (imageA->GetX() != imageB->GetX() ||
      imageA->GetY() != imageB->GetY() || 
      imageA->GetZ() != imageB->GetZ()){
    cerr << "Label images must have equal dimensions" << endl << endl;
    usage();
    exit(1);
  }

  imageA->GetMinMax(&min, &max);
  tmp = max;
  imageB->GetMinMax(&min, &max);
  if (tmp > max)
    max = tmp;

  if (max >= MAX_LABELS){
    cerr << "Maximum number of labels exceeded." << endl << endl;
    usage();
    exit(1);
  }

  //Reset histograms.
  for (i = 0; i < MAX_LABELS; ++i){
    marginalA[i] = 0;
    marginalB[i] = 0;
    for (j = 0; j < MAX_LABELS; ++j){
      jointLabels[i][j] = 0;
    }
  }

  // Reset the list of included labels.
  if (labelFile == NULL){
    for (i = 0; i < MAX_LABELS; ++i){
      includeLabel[i] = 1;
    } 
  } else {
    // First clear all the include labels.
    for (i = 0; i < MAX_LABELS; ++i){
      includeLabel[i] = 0;
    }

    // Read in those to be included.
    ifstream fileIn;
    fileIn.open(labelFile, ifstream::in);
    fileIn.flags(ios_base::skipws);
    while (!fileIn.eof()){
      fileIn >> temp;
      if (temp > 0 && temp < MAX_LABELS){
        includeLabel[temp] = 1;
      }
    }
    fileIn.close();
  }


  pA = imageA->GetPointerToVoxels();
  pB = imageB->GetPointerToVoxels();
  nVoxels = imageA->GetNumberOfVoxels();

  // Fill in histogram.
  for (i = 0; i < nVoxels; ++i){
    if (*pA > 0 || *pB > 0){
      ++marginalA[*pA];
      ++marginalB[*pB];
      ++jointLabels[*pA][*pB];
    }
    ++pA;
    ++pB;
  }


  if (printSummary == true){
    count = 0;
    sumOverlap = sumDice = 0;

    for (i = 1; i <= max; ++i){
      if (includeLabel[i] && (marginalA[i] + marginalB[i] > 0)){
        overlap = jointLabels[i][i] / ((double) ( marginalA[i] + marginalB[i] - jointLabels[i][i] ));
        dice = 2 *  jointLabels[i][i] / ((double) marginalA[i] + marginalB[i] );
        sumOverlap += overlap;
        sumDice += dice;
        ++count;
      }
    }

    if (quiet){
      if (count > 0){
        cout << sumOverlap / ((double) count) << ",";
        cout << sumDice / ((double) count) << endl;
      } else {
        cout << "0,0" << endl;
      }
    } else {
      if (count > 0){
        cout << "Mean OR   : " << sumOverlap / ((double) count) << endl;
        cout << "Mean Dice : " << sumDice / ((double) count) << endl;
      } else {
        cout << "Zero overlap between images." << endl;
      }
    }

  } else {
    // All labels.
    for (i = 1; i <= max; ++i){
      if (includeLabel[i] && (marginalA[i] + marginalB[i] > 0)){

        overlap = jointLabels[i][i] / ((double) ( marginalA[i] + marginalB[i] - jointLabels[i][i] ));
        dice = 2 *  jointLabels[i][i] / ((double) marginalA[i] + marginalB[i] );

        if (diceRow == true){
          cout << dice << ",";
        } else if (orRow == true){
          cout << overlap << ",";
        } else {
          cout << i << "," << marginalA[i] << "," << marginalB[i] << "," << jointLabels[i][i];
          if (printFalseValues == true){
            cout << "," << marginalA[i] - jointLabels[i][i] << "," << marginalB[i] - jointLabels[i][i];
          }
          cout << "," << overlap << "," << dice << endl;
        }

      }
    }

    if (diceRow == true || orRow == true){
      cout << endl;
    }
  }

}
