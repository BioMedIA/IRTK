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

#include <irtkEuclideanDistanceTransform.h>

template <class VoxelType> irtkEuclideanDistanceTransform<VoxelType>::irtkEuclideanDistanceTransform(irtkDistanceTransformMode distanceTransformMode) : irtkImageToImage<VoxelType>()
{
  _distanceTransformMode = distanceTransformMode;
}

/*
 * This procedure computes the squared EDT of a 2D binary image with isotropic
 * voxels of unit dimension. The EDT can obviously be obtained simply from the
 * output by taking the square root of each element. The output can be scaled
 * to account for non-unit dimension. But neither of these functions are
 * provided as options in this procedure.
 *
 * The EDT is returned in the array edt. Memory for the array edt must be
 * allocated by the caller. The binary image can be provided in two different
 * ways. It can be provided in the array img. In this case it is not changed by
 * the procedure. Alternatively, the binary image can be provided in the array
 * edt. In this case, the array img must be the NULL pointer, and the binary
 * image will be overwritten by the EDT. The binary image doesn't have to
 * consist of 0's and 1's. A voxel value of 0 denotes a background voxel and
 * any other value denotes a foreground or feature voxel. The binary image and
 * EDT dimensions are nX x nY voxels.
 *
 * The procedure uses the algorithm described in the paper:
 * CR Maurer Jr, R Qi, V Raghavan. A linear time algorithm for computing exact
 * Euclidean distance transforms of binary images in arbitrary dimensions.
 * IEEE Transactions on Pattern Analysis and Machine Intelligence. In review.
 * A preliminary version of this paper was published in the conference
 * proceedings:
 * CR Maurer Jr, V Raghavan, R Qi. A linear time algorithm for computing the
 * Euclidean distance transform in arbitrary dimensions. In: MF Insana, RM
 * Leahy, eds. Information Processing in Medical Imaging (IPMI) 2001. Berlin:
 * Springer-Verlag, 2001, pp. 358-364. (Davis, CA, June 18-22, 2001).
 */

template <class VoxelType> void irtkEuclideanDistanceTransform<VoxelType>::edtComputeEDT_2D(char *img, long *edt, long nX, long nY)
{
  char *c;
  long i, j, nXY, d, *p, *q, *f;

  /* nXY is number of voxels in 2D image */
  nXY = nX * nY;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != NULL) {
    c = img;
    p = edt;
    for (i = 0; i < nXY; i++, c++, p++) {
      *p = *c;
    }
  }

  /* compute D_1 as simple forward-and-reverse distance propagation */
  /* (instead of calling edtVornoiEDT) */
  /* D_1 is distance to closest feature voxel in row (x direction) */
  /* it is possible to use a simple distance propagation for D_1  because */
  /* L_1 and L_2 norms are equivalent for 1D case */
  for (j = 0; j < nY; j++) {
    /* forward pass */
    p = edt + j * nX;
    d = EDT_MAX_DISTANCE_SQUARED;
    for (i = 0; i < nX; i++, p++) {
      /* set d = 0 when we encounter a feature voxel */
      if (*p) {
        *p = d = 0;
      }
      /* increment distance ... */
      else if (d != EDT_MAX_DISTANCE_SQUARED) {
        *p = ++d;
      }
      /* ... unless we haven't encountered a feature voxel yet */
      else {
        *p = EDT_MAX_DISTANCE_SQUARED;
      }
    }
    /* reverse pass */
    if (*(--p) != EDT_MAX_DISTANCE_SQUARED) {
      d = EDT_MAX_DISTANCE_SQUARED;
      for (i = nX - 1; i >= 0; i--, p--) {
        /* set d = 0 when we encounter a feature voxel */
        if (*p == 0) {
          d = 0;
        }
        /* increment distance after encountering a feature voxel */
        else if (d != EDT_MAX_DISTANCE_SQUARED) {
          /* compare forward and reverse distances */
          if (++d < *p) {
            *p = d;
          }
        }
        /* square distance */
        /* (we use squared distance in rest of algorithm) */
        *p *= *p;
      }
    }
  }

  /* compute D_2 = squared EDT */
  /* solve 1D problem for each column (y direction) */
  f = (long *)malloc(nY * sizeof(long));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_2D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nX; i++) {
    /* fill array f with D_1 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (j = 0; j < nY; j++, p += nX, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT(f, nY)) {
      p = edt + i;
      q = f;
      for (j = 0; j < nY; j++, p += nX, q++) {
        *p = *q;
      }
    }
  }
  free(f);
} /* edtComputeEDT_2D */

template <class VoxelType> void irtkEuclideanDistanceTransform<VoxelType>::edtComputeEDT_3D(char *img, long *edt, long nX, long nY, long nZ)
/*
 * This procedure computes the squared EDT of a 3D binary image with isotropic
 * voxels of unit dimension. See notes for edtComputeEDT_2D.
 */
{
  char *c;
  long i, k, nXY, nXYZ, *p, *q, *f;

  /* nXY is number of voxels in each plane (xy) */
  /* nXYZ is number of voxels in 3D image */
  nXY = nX * nY;
  nXYZ = nX * nY * nZ;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != NULL) {
    c = img;
    p = edt;
    for (i = 0; i < nXYZ; i++, c++, p++) {
      *p = *c;
    }
  }

  /* compute D_2 */
  /* call edtComputeEDT_2D for each plane */
  p = edt;
  for (k = 0; k < nZ; k++, p += nXY) {
    edtComputeEDT_2D(NULL, p, nX, nY);
  }

  /* compute D_3 */
  /* solve 1D problem for each column (z direction) */
  f = (long *)malloc(nZ * sizeof(long));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_3D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nXY; i++) {
    /* fill array f with D_2 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (k = 0; k < nZ; k++, p += nXY, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT(f, nZ)) {
      p = edt + i;
      q = f;
      for (k = 0; k < nZ; k++, p += nXY, q++) {
        *p = *q;
      }
    }
  }
  free(f);
} /* edtComputeEDT_3D */

template <class VoxelType> int irtkEuclideanDistanceTransform<VoxelType>::edtVornoiEDT(long *f, long n)
/*
 * This is Procedure edtVornoiEDT() in tPAMI paper.
 */
{
  long i, l, a, b, c, v, n_S, lhs, rhs;
  static int firstCall = 1;
  static long *g, *h, size_gh = 0;

  /* this procedure is called often */
  /* we keep static arrays to avoid frequent calls to malloc() */
  /* the downside is that this memory is never freed, */
  /* but it is not so much memory: max(nX,nY,nZ)*sizeof(long) */
  /* need to check if arrays are sufficiently large */
  if (n > size_gh && !firstCall) {
    free(g);
    free(h);
    firstCall = 1;
  }
  /* malloc arrays if this is first call to procedure, or if arrays */
  /* are too small and need to be reallocated */
  if (firstCall) {
    g = (long *)malloc(n * sizeof(long));
    if (g == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc g\n");
      exit(EXIT_FAILURE);
    }
    h = (long *)malloc(n * sizeof(long));
    if (h == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc h\n");
      exit(EXIT_FAILURE);
    }
    size_gh = n;
    firstCall = 0;
  }

  /* construct partial Vornoi diagram */
  /* this loop is lines 1-14 in Procedure edtVornoiEDT() in tPAMI paper */
  /* note we use 0 indexing in this program whereas paper uses 1 indexing */
  for (i = 0, l = -1; i < n; i++) {
    /* line 4 */
    if (f[i] != EDT_MAX_DISTANCE_SQUARED) {
      /* line 5 */
      if (l < 1) {
        /* line 6 */
        g[++l] = f[i];
        h[l] = i;
      }
      /* line 7 */
      else {
        /* line 8 */
        while (l >= 1) {
          /* compute removeEDT() in line 8 */
          v = h[l];
          a = v - h[l-1];
          b = i - v;
          c = a + b;
          /* compute Eq. 2 */
          if ((c*g[l] - b*g[l-1] - a*f[i] - a*b*c) > 0) {
            /* line 9 */
            l--;
          } else {
            break;
          }
        }
        /* line 11 */
        g[++l] = f[i];
        h[l] = i;
      }
    }
  }
  /* query partial Vornoi diagram */
  /* this is lines 15-25 in Procedure edtVornoiEDT() in tPAMI paper */
  /* lines 15-17 */
  if ((n_S = l + 1) == 0) {
    return (0);
  }
  /* lines 18-19 */
  for (i = 0, l = 0; i < n; i++) {
    /* line 20 */
    /* we reduce number of arithmetic operations by taking advantage of */
    /* similarities in successive computations instead of treating them as */
    /* independent ones */
    a = h[l] - i;
    lhs = g[l] + a * a;
    while (l < n_S - 1) {
      a = h[l+1] - i;
      rhs = g[l+1] + a * a;
      if (lhs > rhs) {
        /* line 21 */
        l++;
        lhs = rhs;
      } else {
        break;
      }
    }
    /* line 23 */
    /* we put distance into the 1D array that was passed; */
    /* must copy into EDT in calling procedure */
    f[i] = lhs;
  }
  /* line 25 */
  /* return 1 if we queried diagram, 0 if we returned because n_S = 0 */
  return (1);
} /* edtVornoiEDT */


template <class VoxelType> int irtkEuclideanDistanceTransform<VoxelType>::edtVornoiEDT_anisotropic(float *f, long n, float w)
/*
 * This is Procedure edtVornoiEDT() in tPAMI paper.
 */
{
  long i, l, n_S;
  float a, b, c, v, lhs, rhs;
  static int firstCall = 1;
  static long size_gh = 0;
  static float *g, *h;

  /* this procedure is called often */
  /* we keep static arrays to avoid frequent calls to malloc() */
  /* the downside is that this memory is never freed, */
  /* but it is not so much memory: max(nX,nY,nZ)*sizeof(long) */
  /* need to check if arrays are sufficiently large */
  if (n > size_gh && !firstCall) {
    free(g);
    free(h);
    firstCall = 1;
  }
  /* malloc arrays if this is first call to procedure, or if arrays */
  /* are too small and need to be reallocated */
  if (firstCall) {
    g = (float *)malloc(n * sizeof(float));
    if (g == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc g\n");
      exit(EXIT_FAILURE);
    }
    h = (float *)malloc(n * sizeof(float));
    if (h == NULL) {
      fprintf(stderr, "Error in edtVornoiEDT()\n");
      fprintf(stderr, "Cannot malloc h\n");
      exit(EXIT_FAILURE);
    }
    size_gh = n;
    firstCall = 0;
  }

  /* construct partial Vornoi diagram */
  /* this loop is lines 1-14 in Procedure edtVornoiEDT() in tPAMI paper */
  /* note we use 0 indexing in this program whereas paper uses 1 indexing */
  for (i = 0, l = -1; i < n; i++) {
    /* line 4 */
    if (f[i] != EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC) {
      /* line 5 */
      if (l < 1) {
        /* line 6 */
        g[++l] = f[i];
        h[l] = w * i;
      }
      /* line 7 */
      else {
        /* line 8 */
        while (l >= 1) {
          /* compute removeEDT() in line 8 */
          v = h[l];
          a = v - h[l-1];
          b = w * i - v;
          c = a + b;
          /* compute Eq. 2 */
          if ((c*g[l] - b*g[l-1] - a*f[i] - a*b*c) > 0) {
            /* line 9 */
            l--;
          } else {
            break;
          }
        }
        /* line 11 */
        g[++l] = f[i];
        h[l] = w * i;
      }
    }
  }
  /* query partial Vornoi diagram */
  /* this is lines 15-25 in Procedure edtVornoiEDT() in tPAMI paper */
  /* lines 15-17 */
  if ((n_S = l + 1) == 0) {
    return (0);
  }
  /* lines 18-19 */
  for (i = 0, l = 0; i < n; i++) {
    /* line 20 */
    /* we reduce number of arithmetic operations by taking advantage of */
    /* similarities in successive computations instead of treating them as */
    /* independent ones */
    a = h[l] - w * i;
    lhs = g[l] + a * a;
    while (l < n_S - 1) {
      a = h[l+1] - w * i;
      rhs = g[l+1] + a * a;
      if (lhs > rhs) {
        /* line 21 */
        l++;
        lhs = rhs;
      } else {
        break;
      }
    }
    /* line 23 */
    /* we put distance into the 1D array that was passed; */
    /* must copy into EDT in calling procedure */
    f[i] = lhs;
  }
  /* line 25 */
  /* return 1 if we queried diagram, 0 if we returned because n_S = 0 */
  return (1);
} /* edtVornoiEDT_anisotropic */

template <class VoxelType> void irtkEuclideanDistanceTransform<VoxelType>::edtComputeEDT_2D_anisotropic(float *img, float *edt, long nX, long nY, float wX, float wY)
/*
 * This procedure computes the squared EDT of a 2D binary image with anisotropic
 * voxels. See notes for edtComputeEDT_2D. The difference relative to edtComputeEDT_2D
 * is that the edt is a float array instead of a long array, and there are
 * additional parameters for the image voxel dimensions wX and wY.
 */
{
  float *c;
  long i, j, nXY;
  float d, *p, *q, *f;

  /* nXY is number of voxels in 2D image */
  nXY = nX * nY;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != NULL) {
    c = img;
    p = edt;
    for (i = 0; i < nXY; i++, c++, p++) {
      *p = *c;
    }
  }

  /* compute D_1 as simple forward-and-reverse distance propagation */
  /* (instead of calling edtVornoiEDT) */
  /* D_1 is distance to closest feature voxel in row (x direction) */
  /* it is possible to use a simple distance propagation for D_1  because */
  /* L_1 and L_2 norms are equivalent for 1D case */
  for (j = 0; j < nY; j++) {
    /* forward pass */
    p = edt + j * nX;
    d = EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC;
    for (i = 0; i < nX; i++, p++) {
      /* set d = 0 when we encounter a feature voxel */
      if (*p) {
        *p = d = 0;
      }
      /* increment distance ... */
      else if (d != EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC) {
        *p = ++d;
      }
      /* ... unless we haven't encountered a feature voxel yet */
      else {
        *p = EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC;
      }
    }
    /* reverse pass */
    if (*(--p) != EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC) {
      d = EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC;
      for (i = nX - 1; i >= 0; i--, p--) {
        /* set d = 0 when we encounter a feature voxel */
        if (*p == 0) {
          d = 0;
        }
        /* increment distance after encountering a feature voxel */
        else if (d != EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC) {
          /* compare forward and reverse distances */
          if (++d < *p) {
            *p = d;
          }
        }
        /* square distance */
        /* (we use squared distance in rest of algorithm) */
        *p *= wX;
        *p *= *p;
      }
    }
  }

  /* compute D_2 = squared EDT */
  /* solve 1D problem for each column (y direction) */
  f = (float *)malloc(nY * sizeof(float));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_2D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nX; i++) {
    /* fill array f with D_1 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (j = 0; j < nY; j++, p += nX, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT_anisotropic(f, nY, wY)) {
      p = edt + i;
      q = f;
      for (j = 0; j < nY; j++, p += nX, q++) {
        *p = *q;
      }
    }
  }
  free(f);
} /* edtComputeEDT_2D_anisotropic */


/*
 * This procedure computes the squared EDT of a 3D binary image with anisotropic
 * voxels. See notes for edtComputeEDT_2D_anisotropic.
 */

template <class VoxelType> void irtkEuclideanDistanceTransform<VoxelType>::edtComputeEDT_3D_anisotropic(float *img, float *edt, long nX, long nY, long nZ, float wX, float wY, float wZ)
{
  float *c;
  long i, k, nXY, nXYZ;
  float *p, *q, *f;

  /* nXY is number of voxels in each plane (xy) */
  /* nXYZ is number of voxels in 3D image */
  nXY = nX * nY;
  nXYZ = nX * nY * nZ;

  /* if binary image is provided in the array img, copy it to the arry edt */
  /* this is effectively equivalent to computing D_0 */
  if (img != NULL) {
    c = img;
    p = edt;
    for (i = 0; i < nXYZ; i++, c++, p++) {
      *p = *c;
    }
  }

  /* compute D_2 */
  /* call edtComputeEDT_2D for each plane */
  p = edt;
  for (k = 0; k < nZ; k++, p += nXY) {
    edtComputeEDT_2D_anisotropic(NULL, p, nX, nY, wX, wY);
  }

  /* compute D_3 */
  /* solve 1D problem for each column (z direction) */
  f = (float *)malloc(nZ * sizeof(float));
  if (f == NULL) {
    fprintf(stderr, "Error in edtComputeEDT_3D()\n");
    fprintf(stderr, "Cannot malloc f\n");
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nXY; i++) {
    /* fill array f with D_2 distances in column */
    /* this is essentially line 4 in Procedure VoronoiEDT() in tPAMI paper */
    p = edt + i;
    q = f;
    for (k = 0; k < nZ; k++, p += nXY, q++) {
      *q = *p;
    }
    /* call edtVornoiEDT */
    if (edtVornoiEDT_anisotropic(f, nZ, wZ)) {
      p = edt + i;
      q = f;
      for (k = 0; k < nZ; k++, p += nXY, q++) {
        *p = *q *this-> _input->GetZSize() / this->_input->GetXSize();
      }
    }
  }
  free(f);
} /* edtComputeEDT_3D_anisotropic */

template <class VoxelType> void irtkEuclideanDistanceTransform<VoxelType>::Run()
{
  int nx, ny, nz, nt, z, t;
  double wx, wy, wz;

  // Do the initial set up
  this->Initialize();

  // Calculate image dimensions
  nx = this->_input->GetX();
  ny = this->_input->GetY();
  nz = this->_input->GetZ();
  nt = this->_input->GetT();

  // Calculate voxel size
  this->_input->GetPixelSize(&wx, &wy, &wz);
  for ( t = 0; t < nt; t++){
	  if (this->_distanceTransformMode == irtkEuclideanDistanceTransform::irtkDistanceTransform3D) {
		  // Calculate 3D distance transform
		  edtComputeEDT_3D_anisotropic( this->_input->GetPointerToVoxels(0,0,0,t),
			  this->_output->GetPointerToVoxels(0,0,0,t),
			  nx, ny, nz, wx, wy, wz);
	  } else {
		  for (z = 0; z < nz; z++) {
			  // Calculate 2D distance transform slice by slice
			  edtComputeEDT_2D_anisotropic( this->_input->GetPointerToVoxels(0, 0, z, t),
				  this->_output->GetPointerToVoxels(0, 0, z, t),
				  nx, ny, wx, wy);
		  }
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

template <class VoxelType> void irtkEuclideanDistanceTransform<VoxelType>::Radial()
{
  int nx, ny, nz, x, y, z;
  double min;

  // Calculate image dimensions
  nx = this->_input->GetX();
  ny = this->_input->GetY();
  nz = this->_input->GetZ();

  if (this->_distanceTransformMode == irtkEuclideanDistanceTransform::irtkDistanceTransform3D) {
    // Calculate 3D Radial transform
	  min = EDT_MAX_DISTANCE_SQUARED;
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				if(this->_input->GetAsDouble(x,y,z) < min)
					min = this->_input->GetAsDouble(x,y,z);
			  }
		  }
	  }
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				  this->_output->PutAsDouble(x,y,z,this->_input->GetAsDouble(x,y,z)-min);
			  }
		  }
	  }
  } else {
	  // Calculate 2D Radial transform	  
	  for(z=0;z<nz;z++){
		  min = EDT_MAX_DISTANCE_SQUARED;
		  for(y=0;y<ny;y++){
			  for(x=0;x<nz;x++){
				if(this->_input->GetAsDouble(x,y,z) < min)
					min = this->_input->GetAsDouble(x,y,z);
			  }
		  }
		  for(y=0;y<ny;y++){
			  for(x=0;x<nx;x++){
				  this->_output->PutAsDouble(x,y,z,this->_input->GetAsDouble(x,y,z)-min);
			  }
		  }
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

template <class VoxelType> void irtkEuclideanDistanceTransform<VoxelType>::TRadial()
{
  int nx, ny, nz, x, y, z;
  double min;

  // Calculate image dimensions
  nx = this->_input->GetX();
  ny = this->_input->GetY();
  nz = this->_input->GetZ();

  if (this->_distanceTransformMode == irtkEuclideanDistanceTransform::irtkDistanceTransform3D) {
    // Calculate 3D Radial transform
	  min = EDT_MAX_DISTANCE_SQUARED;
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				if(this->_input->GetAsDouble(x,y,z) < min)
					min = this->_input->GetAsDouble(x,y,z);
			  }
		  }
	  }
	  for(x=0;x<nx;x++){
		  for(y=0;y<ny;y++){
			  for(z=0;z<nz;z++){
				  this->_output->PutAsDouble(x,y,z,this->_input->GetAsDouble(x,y,z)-min
					  + abs(this->_input->GetAsDouble(x,y,z))/2);
			  }
		  }
	  }
  } else {
	  // Calculate 2D Radial transform	  
	  for(z=0;z<nz;z++){
		  min = EDT_MAX_DISTANCE_SQUARED;
		  for(y=0;y<ny;y++){
			  for(x=0;x<nz;x++){
				if(this->_input->GetAsDouble(x,y,z) < min)
					min = this->_input->GetAsDouble(x,y,z);
			  }
		  }
		  for(y=0;y<ny;y++){
			  for(x=0;x<nx;x++){
				  this->_output->PutAsDouble(x,y,z,this->_input->GetAsDouble(x,y,z)-min
					  + abs(this->_input->GetAsDouble(x,y,z))/2);
			  }
		  }
	  }
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkEuclideanDistanceTransform<irtkRealPixel>;
