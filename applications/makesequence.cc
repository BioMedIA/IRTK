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

void usage()
{
  cerr << "Usage: makesequence [input 1 ... input n] [output] <options>\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, x, y, z, t;
  irtkImageAttributes ipt0_at;     //add (1/2)
  irtkImageAttributes iptI_at;     //add (1/2)
  
  // Determine how many volumes we have
  t = argc-2;

  if (t < 1) usage();

  cout << "Making sequence from " << t << " volumes" << endl;

  irtkGreyImage* input = new irtkGreyImage[t];

  // Read first image
  cout << "Reading " << argv[1] << endl;
  input[0].Read(argv[1]);
  ipt0_at = input[0].GetImageAttributes();                                                                        //add (1/2)
  
  // Read remaining images
  for (i = 1; i < t; i++) {

    cout << "Reading " << argv[i+1] << endl;
    input[i].Read(argv[i+1]);

  
    //if (!(input[0].GetImageAttributes() == input[i].GetImageAttributes())) {  //removed (1/2)
    //  cerr << "Mismatch of volume geometry" << endl;                          //removed (1/2)
    //  exit(1);                                                                //removed (1/2) -> does not work
    //}                                                                         //removed (1/2)
    
    iptI_at = input[i].GetImageAttributes();                                                                        //add (1/2)
    if ((ipt0_at._x!=iptI_at._x)||(ipt0_at._y!=iptI_at._y)||(ipt0_at._z!=iptI_at._z)||(ipt0_at._t!=iptI_at._t)){    //add (1/2)
	    cerr << "Mismatch of volume geometry" << endl;                                                          //add (1/2)
	    exit(1);                                                                                                //add (1/2)
    }                                                                                                               //add (1/2)
  }
  
  //irtkGreyImage output(input[0].GetImageAttributes());              //removed (2/2)
  
  irtkImageAttributes OutputSeqAttributes=input[0].GetImageAttributes();  //add (2/2)
  OutputSeqAttributes._t=t;                                               //add (2/2)
  irtkGreyImage output(OutputSeqAttributes);                              //add (2/2)
  

  cout << "Inserting volumes into sequence" << endl;
  for (i = 0; i < t; i++) {
    cout << "Volume " << i+1 << " ..." << endl;
    for (z = 0; z < output.GetZ(); z++) {
      for (y = 0; y < output.GetY(); y++) {
        for (x = 0; x < output.GetX(); x++) {
          output(x, y, z, i) = input[i](x, y, z);
        }
      }
    }
  }

  // Write image
  cout << "Writing sequence to " << argv[t+1] << endl;
  output.Write(argv[t+1]);

  delete[] input;
}

