#ifndef _EXTERNAL_LIBS_H
#define _EXTERNAL_LIBS_H

#include "common.h"

//------------------------------------ MPI functions ---------------------------------------------
extern "C" {
   void gsyncX();
   void gabcastX(int type, char* x, int n, int root);
   void gdbcastX(int type, double* x, int n, int root);
   void gibcastX(int type, int* x, int n, int root);
   void gdsumX(int type, double* x, int n);
};
//------------------------------------------------------------------------------------------------


//------------------------------- fortran functions ----------------------------------------------
extern "C" {

   double truesdvel1_(int* model, double *x, double *y, double *z);
   double truesdvel2_(int* model, double *x, double *y, double *z);
   double sdparameters_(int *pr);

   void getnbrsize_(int *mynbr, int *mynbrSize, int &n1, int &n2, int &n3);

   void copyiface_(int *mynbr, int &n1, int &n2, int &n3,
                   double *bf1, double *bf2, double *bf3,
                   double *bf1_copy, double *bf2_copy, double *bf3_copy);

   void swapface_(int *mynbr, int *mynbrSize, int &n1, int &n2, int &n3,
		  double *bf1, double *bf2, double *bf3, double *df1, double *df2, double *df3,
		  int &s, double *a1, double *a2, int &FLUXTYPE, int &BdryVarType,
		  double*, double*, double*, int &mortarflag, int *mortartype, int*,
		  double*, double*, double*, double*, double*, double*,
		  double*, double*, double*, double*, double*, double*,
		  double*, double*, double*, double*, double*, double*,
		  int &AvgType, int*, int*, int*, int*, int*, int *,
		  int*, int*, int*, int*, int*, int*,
		  int*, int*, int*, int*, int*, int*,
		  double*, double*, double*, double*, double*, double*);



  // LAPACK
  void dgesv_(int& n, int& nrhs, double* A, int& lda, int* ipiv, double* b, int& ldb, int& info);
  void dposv_(char& uplo, int& n, int& nrhs, double* A, int& lda, double* B, int& ldb, int& info);
};

#endif
