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

#include <irtkGIPL.h>

int main(int argc, char **argv)
{
  unsigned int magic_number;

  // Open file
  ifstream from(argv[1]);

  // Read magic no
  from.seekg(252, ios::beg);
  from.read((char *)&magic_number, 4);
  from.close();

#ifndef WORDS_BIGENDIAN
  swap32((char *)&magic_number, (char *)&magic_number, 1);
#endif

  // Check header
  if (magic_number == GIPL_MAGIC1) {
    cout << "File has Guy's magic no., writing UNC magic no." << endl;
    magic_number = GIPL_MAGIC2;
  } else {
    if (magic_number == GIPL_MAGIC2) {
      cout << "File has UNC magic no., writing Guy's magic no." << endl;
      magic_number = GIPL_MAGIC1;
    } else {
      cout << "File has not UNC or Guy's magic no." << endl;
      exit(1);
    }
  }

  ofstream to(argv[1], ofstream::ate);

#ifndef WORDS_BIGENDIAN
  swap32((char *)&magic_number, (char *)&magic_number, 1);
#endif

  // Write magic no
  to.seekp(252, ios::beg);
  to.write((char *)&magic_number, 4);
  to.close();
}






