%define DOCSTRING
"This module contains wrapper classes for the `Image Registration Toolkit` (ITK) library.

For more information about the ITK library see:

http://www.doc.ic.ac.uk/~dr/software/"
%enddef

%module(docstring=DOCSTRING) itk
%{
  #include <itkCommon.h>
  #include <itkGeometry.h>
  #include <itkImage.h>
  #include <itkTransformation.h>
  #include <itkLatticeFreeFormTransformation.h>
  #include <itkMultiFrameLatticeFreeFormTransformation.h>
  #include <itkCardiac.h>

  #include <itkImageToImage.h>
  #include <itkGaussianBlurring.h>
  #include <itkVector3D.h>

  #include <iomanip>
%}

%include "typemaps.i"

%typemap(in) double DOUBLE3INPUT[3] (double double3InputCopy[3])
{
  int i;
  if (!PySequence_Check($input))
  {
    PyErr_SetString(PyExc_ValueError, "Expected a sequence.");
    return NULL;
  }
  if (PySequence_Length($input) != 3)
  {
    PyErr_SetString(PyExc_ValueError, "Size mismatch. Expected 3 elements.");
    return NULL;
  }
  for (i = 0; i < 3; i++)
  {
    PyObject* o = PySequence_GetItem($input, i);
    if (PyNumber_Check(o))
    {
      double3InputCopy[i] = (double) PyFloat_AsDouble(o);
    }
    else
    {
      PyErr_SetString(PyExc_ValueError, "Sequence elements must be numbers.");
      return NULL;
    }
  }
  $1 = double3InputCopy;
}

%typemap(argout) double DOUBLE3OUTPUT[3]
{
  PyObject* o;
  PyObject* o2;
  PyObject* o3;

  o = PyList_New(3);
  PyList_SetItem(o, 0, PyFloat_FromDouble($1[0]));
  PyList_SetItem(o, 1, PyFloat_FromDouble($1[1]));
  PyList_SetItem(o, 2, PyFloat_FromDouble($1[2]));

  if ((!$result) || ($result == Py_None))
  {
    $result = o;
  }
  else
  {
    if (!PyTuple_Check($result))
    {
      PyObject *o2 = $result;
      $result = PyTuple_New(1);
      PyTuple_SetItem($result, 0, o2);
    }
    o3 = PyTuple_New(1);
    PyTuple_SetItem(o3, 0, o);
    o2 = $result;
    $result = PySequence_Concat(o2, o3);
    Py_DECREF(o2);
    Py_DECREF(o3);
  }
}

%typemap(in, numinputs=0) double DOUBLE3OUTPUT[3] (double temp[3])
{
  $1 = temp;
}

%typemap(in) (int NINTARGS, int* INTARGS)
{
  int i;
  if (!PyList_Check($input))
  {
    PyErr_SetString(PyExc_ValueError, "Expected a list.");
    return NULL;
  }
  $1 = PyList_Size($input);
  $2 = new int[$1];
  for (i = 0; i < $1; i++)
  {
    PyObject* o = PyList_GetItem($input, i);
    if (!PyInt_Check(o))
    {
      delete[] $2;
      PyErr_SetString(PyExc_ValueError, "List items must be integers.");
      return NULL;
    }
    $2[i] = PyInt_AsInt(o);
  }
}

%typemap(freearg) (int NINTARGS, int* INTARGS)
{
  if ($2)
    delete[] $2;
}

%include "geometry++/geometry.i"
%include "image++/image.i"
%include "packages/transformation++/transformation.i"
%include "packages/cardiac/cardiac.i"
