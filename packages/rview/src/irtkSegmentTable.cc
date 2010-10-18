/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkSegmentTable.h>

irtkSegmentTable::irtkSegmentTable()
{
}

irtkSegmentTable::~irtkSegmentTable()
{
}

void irtkSegmentTable::Set(int id, char* label, unsigned char r, unsigned char g, unsigned char b, double trans, int vis)
{
  _entry[id].setLabel(label);
  _entry[id].setColor(r, g, b);
  _entry[id].setTrans(trans);
  _entry[id].setVisibility(vis);
}

void irtkSegmentTable::SetLabel(int id, char* label)
{
  _entry[id].setLabel(label);
}

void irtkSegmentTable::SetColor(int id, unsigned char red, unsigned char green, unsigned char blue)
{
  _entry[id].setColor(red, green, blue);
}

void irtkSegmentTable::SetTrans(int id, double t)
{
  _entry[id].setTrans(t);
}

void irtkSegmentTable::SetVisibility(int id, int vis)
{
  _entry[id].setVisibility(vis);
}

char *irtkSegmentTable::Get(int id, unsigned char* r, unsigned char* g, unsigned char* b, double* trans, int* v) const
{
  // Get r,g,b
  _entry[id].getColor(r, g, b);

  // Get transparency
  *trans = _entry[id].getTrans();

  // Get visibility
  *v = _entry[id].getVisibility();

  return _entry[id].getLabel();
}

void irtkSegmentTable::Clear()
{
  int i;

  for (i = 0; i < this->Size(); i++) {
    if (IsValid(i) == true) Clear(i);
  }
}

void irtkSegmentTable::Clear(int id)
{
  _entry[id].setLabel(NULL);
  _entry[id].setColor(0, 0, 0);
  _entry[id].setTrans(0);
}

void irtkSegmentTable::Read(char *name)
{
  double trans;
  char c, buffer[256];
  int i, n, r, g, b, id, vis;

  // Clear entries
  Clear();

  ifstream from(name);
  if (!from) {
    cerr << "irtkSegmentTable::Read: Can't open file " << name << "\n";
    exit(1);
  }

  // Read keyword
  from >> buffer;
  if ((strcmp(buffer, "irtkSegmentTable:") != 0) && (strcmp(buffer, "itkSegmentTable:") != 0)) {
    cerr << "irtkSegmentTable::Read: Not a valid segment table" << endl;
    exit(1);
  }

  // Read number of valid entries
  from >> n;

  // Read entries
  for (i = 0; i < n; i++) {
    from >> id >> r >> g >> b >> trans >> vis;
    for (;;) {
      c = from.peek();
      if ((c == ' ') || (c == '\t')) {
        from.get(c);
      } else {
        break;
      }
    }
    from.getline(buffer, 255);
    _entry[id].setLabel(buffer);
    _entry[id].setColor(r, g, b);
    _entry[id].setTrans(trans);
    _entry[id].setVisibility(vis);
  }
}

void irtkSegmentTable::Write(char *name)
{
  int id, n;
  unsigned char r, g, b;

  // Open file
  ofstream to(name);

  if (!to) {
    cerr << "irtkSegmentTable::Write: Can't open file " << name << "\n";
    exit(1);
  }

  // Count number of valid entries
  n = 0;
  for (id = 0; id < this->Size(); id++) {
    if (IsValid(id) == true) n++;
  }

  // Write header
  to << "irtkSegmentTable: " << n << endl;

  // Write entries
  for (id = 0; id < this->Size(); id++) {
    if (IsValid(id) == true) {
      _entry[id].getColor(&r, &g, &b);
      to << id << "\t" << int(r) << "\t" << int(g) << "\t" << int(b) << "\t" << _entry[id].getTrans()<< "\t" << _entry[id].getVisibility() << "\t" << _entry[id].getLabel() << endl;
    }
  }
}


