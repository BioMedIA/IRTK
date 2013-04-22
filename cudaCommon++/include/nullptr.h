/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    :
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      :
  Version   :
  Changes   : $Author: bkainz $
  Original reduction implementation by: Noel Lopes, GPUMLib

=========================================================================*/


const                        // this is a const object...
class {
public:
  template<class T>          // convertible to any type
    operator T*() const      // of null non-member
    { return 0; }            // pointer...
  template<class C, class T> // or any type of null
    operator T C::*() const  // member pointer...
    { return 0; }
private:
  void operator&() const;    // whose address can't be taken
} nullptr = {};