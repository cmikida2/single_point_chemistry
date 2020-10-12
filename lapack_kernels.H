#ifndef __LAPACK_KERNELS_H__
#define __LAPACK_KERNELS_H__

//#include "FC.h"

extern "C" {

  void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
  void dgesvd_(char*, char*, int*, int*, double [], int*, double*, double*, int*, double*, int*, double [], int*, int*);
  double dnrm2_(int *, double *, int *);

}
#endif
