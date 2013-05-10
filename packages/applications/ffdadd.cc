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

#include <irtkTransformation.h>

char *target_name = NULL;
char *dofin_name  = NULL;
char *dofout_name = NULL;

irtkGreyImage *target;
irtkMultiLevelFreeFormTransformation *mffd;

#ifdef WIN32
#include <time.h>
#include <windows.h>

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone 
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
    FILETIME ft;
    unsigned __int64 tmpres = 0;
    static int tzflag;

    if (NULL != tv)
    {
        GetSystemTimeAsFileTime(&ft);

        tmpres |= ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |= ft.dwLowDateTime;

        /*converting file time to unix epoch*/
        tmpres -= DELTA_EPOCH_IN_MICROSECS; 
        tmpres /= 10;  /*convert into microseconds*/
        tv->tv_sec = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }

    if (NULL != tz)
    {
        if (!tzflag)
        {
            _tzset();
            tzflag++;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime = _daylight;
    }

    return 0;
}
#else
#include <sys/time.h>
#endif

#include <nr.h>
#include <nrutil.h>

void usage()
{
  cerr << "Usage: ffdadd [target] [dofin] [dofout] [options]" << endl;
  cerr << "where the following options are supported:\n" << endl;
  cerr << "<-Tx1 value>  Region of interest in target image" << endl;
  cerr << "<-Ty1 value>  Region of interest in target image" << endl;
  cerr << "<-Tz1 value>  Region of interest in target image" << endl;
  cerr << "<-Tt1 value>  Region of interest in target image" << endl;
  cerr << "<-Tx2 value>  Region of interest in target image" << endl;
  cerr << "<-Ty2 value>  Region of interest in target image" << endl;
  cerr << "<-Tz2 value>  Region of interest in target image" << endl;
  cerr << "<-Tt2 value>  Region of interest in target image" << endl;
  cerr << "<-dx  value>  Control point spacing for x" << endl;
  cerr << "<-dy  value>  Control point spacing for y" << endl;
  cerr << "<-dz  value>  Control point spacing for z" << endl;
  cerr << "<-dt  value>  Control point spacing for t" << endl;
  cerr << "<-ds  value>  Control point spacing for x, y, z and t" << endl;
  cerr << "<-3D>         Generate 3D FFD (default)" << endl;
  cerr << "<-4D>         Generate 4D FFD" << endl;
  cout << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double dx, dy, dz, dt, noise;
  int i1, j1, k1, l1, i2, j2, k2, l2, ok, fluid, ffd4D;

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

  // Read image image
  cout << "Reading image ... "; cout.flush();
  target = new irtkGreyImage;
  target->Read(target_name);
  cout << "done" << endl;

  // Fix ROI
  i1 = 0;
  j1 = 0;
  k1 = 0;
  l1 = 0;
  i2 = target->GetX();
  j2 = target->GetY();
  k2 = target->GetZ();
  l2 = target->GetT();

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;
  dt = 20;

  // No 4D FFD
  ffd4D = false;

  // No fluid
  fluid = false;

  // Don't add noise
  noise = 0;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-Tx1") == 0)) {
      argc--;
      argv++;
      i1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx2") == 0)) {
      argc--;
      argv++;
      i2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty1") == 0)) {
      argc--;
      argv++;
      j1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty2") == 0)) {
      argc--;
      argv++;
      j2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz1") == 0)) {
      argc--;
      argv++;
      k1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz2") == 0)) {
      argc--;
      argv++;
      k2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tt1") == 0)) {
      argc--;
      argv++;
      l1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tt2") == 0)) {
      argc--;
      argv++;
      l2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dx") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dy") == 0)) {
      argc--;
      argv++;
      dy = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dz") == 0)) {
      argc--;
      argv++;
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dt") == 0)) {
      argc--;
      argv++;
      dt = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      dt = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-fluid") == 0)) {
      argc--;
      argv++;
      fluid = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-3D") == 0)) {
      argc--;
      argv++;
      ffd4D = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-4D") == 0)) {
      argc--;
      argv++;
      ffd4D = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-noise") == 0)) {
      argc--;
      argv++;
      noise = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read transformation
  if (fluid == true) {
    mffd = new irtkFluidFreeFormTransformation;
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }
  mffd->irtkTransformation::Read(dofin_name);

  // If there is an region of interest, use it
  if ((i1 != 0) || (i2 != target->GetX()) ||
      (j1 != 0) || (j2 != target->GetY()) ||
      (k1 != 0) || (k2 != target->GetZ()) ||
      (l1 != 0) || (l2 != target->GetT())) {
    *target = target->GetRegion(i1, j1, k1, l1, i2, j2, k2, l2);
  }

  if (ffd4D == false) {
    // Create free-form transformation
    irtkBSplineFreeFormTransformation3D *affd = new irtkBSplineFreeFormTransformation(*target, dx, dy, dz);

    // Add noise to the transformation
    if (noise > 0) {
      int i;
      long temp;

      timeval tv;
      gettimeofday(&tv, NULL);
      temp = -tv.tv_usec;

      for (i = 0; i < affd->NumberOfDOFs(); i++) affd->Put(i, noise*gasdev(&temp));
    }

    // Add and write file
    mffd->PushLocalTransformation(affd);
  } else {
    // Create free-form transformation
    irtkBSplineFreeFormTransformation4D *affd = new irtkBSplineFreeFormTransformation4D(*target, dx, dy, dz, dt);

    // Add noise to the transformation
    if (noise > 0) {
      int i;
      long temp;

      timeval tv;
      gettimeofday(&tv, NULL);
      temp = -tv.tv_usec;

      for (i = 0; i < affd->NumberOfDOFs(); i++) affd->Put(i, noise*gasdev(&temp));
    }

    // Add and write file
    mffd->PushLocalTransformation(affd);
  }

  mffd->irtkTransformation::Write(dofout_name);
}
