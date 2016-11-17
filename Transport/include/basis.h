#ifndef _BASIS_H
#define _BASIS_H

#include "common.h"

//------------------------------------- local basis ----------------------------------------------

/*

 local ordering of the nodes:

  Y                      3______2
                        /|      /
  |                    / |     /|
  |                   /7_|___6/ |
  |                   |  |    | |
  |______ X           |  0____|_1
  /                   | /     | /
 /                    |/      |/
Z                     /4_____5/

*/
// basis functions on a reference cube [-1,1] x [-1,1] x [-1,1]

inline double w0(double x, double y, double z) {return ( 1.0 - x ) * ( 1.0 - y ) * ( 1.0 - z ) / 8.0;}
inline double w1(double x, double y, double z) {return ( 1.0 + x ) * ( 1.0 - y ) * ( 1.0 - z ) / 8.0;}
inline double w2(double x, double y, double z) {return ( 1.0 + x ) * ( 1.0 + y ) * ( 1.0 - z ) / 8.0;}
inline double w3(double x, double y, double z) {return ( 1.0 - x ) * ( 1.0 + y ) * ( 1.0 - z ) / 8.0;}
inline double w4(double x, double y, double z) {return ( 1.0 - x ) * ( 1.0 - y ) * ( 1.0 + z ) / 8.0;}
inline double w5(double x, double y, double z) {return ( 1.0 + x ) * ( 1.0 - y ) * ( 1.0 + z ) / 8.0;}
inline double w6(double x, double y, double z) {return ( 1.0 + x ) * ( 1.0 + y ) * ( 1.0 + z ) / 8.0;}
inline double w7(double x, double y, double z) {return ( 1.0 - x ) * ( 1.0 + y ) * ( 1.0 + z ) / 8.0;}


// partial derivatives of the basis functions

inline double w0x(double x, double y, double z) {return  -  ( 1.0 - y ) * ( 1.0 - z ) / 8.0;}
inline double w0y(double x, double y, double z) {return  -  ( 1.0 - x ) * ( 1.0 - z ) / 8.0;}
inline double w0z(double x, double y, double z) {return  -  ( 1.0 - x ) * ( 1.0 - y ) / 8.0;}

inline double w1x(double x, double y, double z) {return     ( 1.0 - y ) * ( 1.0 - z ) / 8.0;}
inline double w1y(double x, double y, double z) {return  -  ( 1.0 + x ) * ( 1.0 - z ) / 8.0;}
inline double w1z(double x, double y, double z) {return  -  ( 1.0 + x ) * ( 1.0 - y ) / 8.0;}

inline double w2x(double x, double y, double z) {return     ( 1.0 + y ) * ( 1.0 - z ) / 8.0;}
inline double w2y(double x, double y, double z) {return     ( 1.0 + x ) * ( 1.0 - z ) / 8.0;}
inline double w2z(double x, double y, double z) {return  -  ( 1.0 + x ) * ( 1.0 + y ) / 8.0;}

inline double w3x(double x, double y, double z) {return  -  ( 1.0 + y ) * ( 1.0 - z ) / 8.0;}
inline double w3y(double x, double y, double z) {return     ( 1.0 - x ) * ( 1.0 - z ) / 8.0;}
inline double w3z(double x, double y, double z) {return  -  ( 1.0 - x ) * ( 1.0 + y ) / 8.0;}

inline double w4x(double x, double y, double z) {return  -  ( 1.0 - y ) * ( 1.0 + z ) / 8.0;}
inline double w4y(double x, double y, double z) {return  -  ( 1.0 - x ) * ( 1.0 + z ) / 8.0;}
inline double w4z(double x, double y, double z) {return     ( 1.0 - x ) * ( 1.0 - y ) / 8.0;}

inline double w5x(double x, double y, double z) {return     ( 1.0 - y ) * ( 1.0 + z ) / 8.0;}
inline double w5y(double x, double y, double z) {return  -  ( 1.0 + x ) * ( 1.0 + z ) / 8.0;}
inline double w5z(double x, double y, double z) {return     ( 1.0 + x ) * ( 1.0 - y ) / 8.0;}

inline double w6x(double x, double y, double z) {return     ( 1.0 + y ) * ( 1.0 + z ) / 8.0;}
inline double w6y(double x, double y, double z) {return     ( 1.0 + x ) * ( 1.0 + z ) / 8.0;}
inline double w6z(double x, double y, double z) {return     ( 1.0 + x ) * ( 1.0 + y ) / 8.0;}

inline double w7x(double x, double y, double z) {return  -  ( 1.0 + y ) * ( 1.0 + z ) / 8.0;}
inline double w7y(double x, double y, double z) {return     ( 1.0 - x ) * ( 1.0 + z ) / 8.0;}
inline double w7z(double x, double y, double z) {return     ( 1.0 - x ) * ( 1.0 + y ) / 8.0;}


// define a type for a basis function pointer 
typedef double (*function)(double, double, double);

//------------------------------------------------------------------------------------------------

// compute det(JF) where  F  :  Eref   ->   E
inline double det_JF(double x[], double y[], double z[])
{ return ( x[1]-x[0] ) * ( y[1]-y[0] ) * ( z[1]-z[0] ) / 8.0  ; }

inline double det_JF(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{ return (xmax-xmin) * (ymin-ymax) * (zmin-zmax) / 8.0 ; }

//------------------------------------------------------------------------------------------------

// identity function
inline double one(double x, double y, double z) { return 1.0; }

//------------------------------------------------------------------------------------------------

#endif
