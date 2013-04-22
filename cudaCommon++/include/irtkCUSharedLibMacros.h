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


#ifndef irtkCUSharedLibMacros_H_
#define irtkCUSharedLibMacros_H_

//TODO: with IRTK as static base lib, it will be hard to get a shared plugin ready system
//TODO: make IRTK to be built as shared libs

#ifdef WIN32
#pragma warning( disable : 4251 ) // disable the warning about exported template code from stl

   #ifdef irtkCUcommonLib_USE_STATIC
      #define irtkCULib_DLLAPI
   #else
      #ifdef irtkCUcommonLib_EXPORTS
         #define irtkCULib_DLLAPI __declspec(dllexport)
      #else
         #define irtkCULib_DLLAPI __declspec(dllimport)
      #endif
   #endif
#else
   #define irtkCULib_DLLAPI
#endif

#endif // irtkCUSharedLibMacros_H_
