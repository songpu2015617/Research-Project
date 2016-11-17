#ifndef _UTILITIES_H
#define _UTILITIES_H

#include "common.h"
#include "external_libs.h"
#include "eqn.h"
#include "io.h"
#include "basis.h"
#include "quad.h"
#include "problem.h"

//--------------- solve an nxn linear system using Cholesky factorization ------------------------
// A  x    =   U^t U  x   = b
// Note: U is upper triangular
void cholesky_solve(int n, double U[], double x[], double b[]);
//------------------------------------------------------------------------------------------------

//------------------------ compute the determinant of a 3x3 matrix -------------------------------
inline double det(double A[][3]) {
   return  A[0][0] * ( A[1][1]*A[2][2] - A[2][1]*A[1][2] )
         - A[0][1] * ( A[1][0]*A[2][2] - A[2][0]*A[1][2] )
         + A[0][2] * ( A[1][0]*A[2][1] - A[2][0]*A[1][1] );
}
//------------------------------------------------------------------------------------------------

//---------------------- solve a 3x3 linear system using Cramer's rule ---------------------------
void cramer(double A[][3], double b[], double x[]);
//------------------------------------------------------------------------------------------------


//---------------- solve a 3x3 linear system using Cholesky factorization ------------------------
void cholesky_factor(double A[][3]); // not implemented
void cholesky_solve(double A[][3], double b[], double x[]);
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------


//---------------------------- compute the subdomain averages ------------------------------------
void avr_x(double c[][8], double cbarx[][8], int nx, int ny, int nz); // average along x-
void avr_y(double c[][8], double cbary[][8], int nx, int ny, int nz); // average along y-
void avr_z(double c[][8], double cbarz[][8], int nx, int ny, int nz); // average along z-
//------------------------------------------------------------------------------------------------


//------- functions to handle the temporary buffers for the interface values of c and cbar -------
void load_buffers(double c[][8], double bf1[], double bf2[], double bf3[], int nx, int ny, int nz);
void read_buffers(double cbarx[][8], double cbary[][8], double cbarz[][8],
                  double bf1[], double bf2[], double bf3[],
                  int nx, int ny, int nz);
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
void swapbdry( int n1, int n2, int n3, int *mynbr, int mynbrSize[],
               double *bf1, double *bf2, double *bf3, double *df1, double *df2, double *df3,
               int BdryVarType, int AvgType);
//------------------------------------------------------------------------------------------------


//----------------- get the the interface values from the neighbouring processor -----------------
void interface_val( double c[][8], double upwd1[][8], double upwd2[][8], double upwd3[][8],
                    int n1, int n2, int n3, int *mynbr, fstream& out );
//------------------------------------------------------------------------------------------------


//---------------- compute the averaged concentration in the whole domain ------------------------
void domain_avr_scalar(double c[][8], double cbarx[][8], double cbary[][8], double cbarz[][8],
                       int n1, int n2, int n3, int *mynbr, fstream& out);
//------------------------------------------------------------------------------------------------


//-------------------- compute the averaged flux in the whole domain -----------------------------
void domain_avr_vector(double qx[][8], double qy[][8], double qz[][8],
                       double qbarx[][8], double qbary[][8], double qbarz[][8],
                       int n1, int n2, int n3, int *mynbr, fstream& out);
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
void compute_flux( eqn *Ptr, int *mynbr, int &model,
                   double c[][8], double cbarx[][8], double cbary[][8], double cbarz[][8],
                   double qx[][8], double qy[][8], double qz[][8],
                   int& dflag, int& dtype, fstream& out );
//------------------------------------------------------------------------------------------------


//---------------------------- calculate the L2 norm of the error --------------------------------
double conc_L2error_sq_subd(int& model, eqn *Ptr, double c[][8], double &t, int& testnumber);
double conc_L2error_sq(int& model, eqn *Ptr, double c[][8], double &t, int& testnumber);
double conc_L2error(int& model, eqn *Ptr, double c[][8], double &t, int& testnumber);
double flux_L2error_sq_subd(int& model, eqn *Ptr, double q[][8], int& dir, double &t, int& testnumber);
double flux_L2error_sq(int& model, eqn *Ptr, double qx[][8], double qy[][8], double qz[][8],
                       double &t, int& testnumber);
double flux_L2error(int& model, eqn *Ptr, double qx[][8], double qy[][8], double qz[][8],
                    double &t, int& testnumber);
//------------------------------------------------------------------------------------------------


//--------------- compute the mean and the variance of the concentration -------------------------
void compute_stochparam(int n1, int n2, int n3,
                        double c[][8], double &cweight, fstream& stochIO);
bool isRealizationToPlot(int& StochCount, int *RealizationIdx, int *NumOfRealizations);
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
double mass(eqn *Ptr, double c[][8]);
//------------------------------------------------------------------------------------------------

#endif
