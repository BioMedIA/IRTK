/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCommon.h>

#define MAX_LINE 255

int ReadInt(ifstream &in)
{
  char c;
  char *string, *s;
  int data;

  string = new char[MAX_LINE];
  s = string;
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadInt: No valid line format\n";
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  data = atoi(string);
  delete s;
  return (data);
}

//
// General routine to read float from a file stream.
//
float ReadFloat(ifstream &in)
{
  char c;
  char *string, *s;
  float data;

  string = new char[MAX_LINE];
  s = string;
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadFloat: No valid line format\n";
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  data = atof(string);
  delete s;
  return (data);
}

// General routine to read list of char (string) from a file stream.
char *ReadString(ifstream &in)
{
  char c;
  char *string;

  string = new char[MAX_LINE];
  do {
    in.get(string, MAX_LINE, '\n');
    in.get(c);
  } while ((strlen(string) == 0) || (string[0] == '#'));
  if ((string = strchr(string, '=')) == NULL) {
    cerr << "ReadString: No valid line format\n";
    exit(1);
  }
  do {
    string++;
  } while ((*string == ' ') || (*string == '\t'));

  return (string);
}
