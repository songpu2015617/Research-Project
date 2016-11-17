#ifndef _LIMITER_H
#define _LIMITER_H

#include "common.h"
#include "utilities.h"
#include "external_libs.h"

//------------------------------------------------------------------------------------------------
inline int sign(int var) {  if (var!=0) return (var>0) ? 1 : -1; else return 0;  }
inline int sign(double var) {  if( fabs(var)>1.0e-16 ) return (var>0.0) ? 1 : -1; else return 0;  }
double minmod(double a1, double a2, double a3);
inline double arithmetic_mean(double c1, double c2, double c3) { return ( c1 + c2 + c3 ) / 3.0; }
inline double arithmetic_mean(double c1, double c2, double c3, double c4) { return ( c1 + c2 + c3 + c4) / 4.0; }
double cell_avr(double c[8]);
void cell_avr_array(int& n1, int& n2, int& n3, double c[][8], double avr_c[]);
void limiter(double& alpha, int& n1, int& n2, int& n3, int *mynbr, double c[][8], fstream& out);
//------------------------------------------------------------------------------------------------

#endif
