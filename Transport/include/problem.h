#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "common.h"
#include "external_libs.h"
#include "eqn.h"



const int velN = 14;  // this how many analytical velocity fields are available

const double K1 = 1.0, K2 = 1.0, x0 = 0.1; // used for test cases 3 and 4

const double ax = 0.5, ay = -0.5;          // used for test case 9

const double Ux = 0.2; // const velocity used for some of the 1D test cases

inline double c_inflow_X(double x, double y, double z) { return ( (fabs(x-0.0) < 1.0e-16) ? 1.0 : 0.0 ); }
inline double c_inflow_Y(double x, double y, double z) { return ( (fabs(y-0.0) < 1.0e-16) ? 1.0 : 0.0 ); }
inline double c_inflow_Z(double x, double y, double z) { return ( (fabs(z-0.0) < 1.0e-16) ? 1.0 : 0.0 ); }
double c_inflow(int& model, double x, double y, double z, double t, int testnumber);

double c_initial(int& model, double x, double y, double z, int testnumber);
double f(int& model, double x, double y, double z, double t, int testnumber, eqn *Ptr);

double c_true(int& model, double x, double y, double z, double t, int testnumber);
double qx_true(int& model, eqn *Ptr, double x, double y, double z, double t, int testnumber);
double qy_true(int& model, eqn *Ptr, double x, double y, double z, double t, int testnumber);
double qz_true(int& model, eqn *Ptr, double x, double y, double z, double t, int testnumber);

//---------------------------------  analytical velocity -----------------------------------------
// the velocity functions are called by
// c_initial(), f(), c_true(), qx_true(), qy_true(), qz_true(), and setvel_()

// define a type for a velocity function pointer 
typedef double (*velfun)(int&, double, double, double);

inline double trueU1_test1(int& model, double x, double y, double z) { return Ux;  } // 1D convergence test
inline double trueU2_test1(int& model, double x, double y, double z) { return 0.0; }
inline double trueU3_test1(int& model, double x, double y, double z) { return 0.0; }

inline double trueU1_test2(int& model, double x, double y, double z) { return Ux;  } // 1D convergence test
inline double trueU2_test2(int& model, double x, double y, double z) { return 0.0; }
inline double trueU3_test2(int& model, double x, double y, double z) { return 0.0; }

inline double trueU1_test3(int& model, double x, double y, double z) { return Ux;  } // 1D convergence test
inline double trueU2_test3(int& model, double x, double y, double z) { return 0.0; }
inline double trueU3_test3(int& model, double x, double y, double z) { return 0.0; }

inline double trueU1_test4(int& model, double x, double y, double z) { return Ux;  } // 1D convergence test
inline double trueU2_test4(int& model, double x, double y, double z) { return 0.0; }
inline double trueU3_test4(int& model, double x, double y, double z) { return 0.0; }

inline double trueU1_test5(int& model, double x, double y, double z) { return Ux;  } // 1D convergence test
inline double trueU2_test5(int& model, double x, double y, double z) { return 0.0; }
inline double trueU3_test5(int& model, double x, double y, double z) { return 0.0; }

inline double trueU1_test6(int& model, double x, double y, double z) { return Ux;  } // 1D convergence test
inline double trueU2_test6(int& model, double x, double y, double z) { return 0.0; }
inline double trueU3_test6(int& model, double x, double y, double z) { return 0.0; }

inline double trueU1_test7(int& model, double x, double y, double z) { return Ux;  } // sharp 1-D front
inline double trueU2_test7(int& model, double x, double y, double z) { return 0.0; } // propagating along X
inline double trueU3_test7(int& model, double x, double y, double z) { return 0.0; } // ( no diffusion )

inline double trueU1_test8(int& model, double x, double y, double z) { return 0.0; } // used to test
inline double trueU2_test8(int& model, double x, double y, double z) { return 0.0; } // the LDG scheme 
inline double trueU3_test8(int& model, double x, double y, double z) { return 0.0; } // with diffusion only

inline double trueU1_test9(int& model, double x, double y, double z) { return ax;  } // these are to test for convergence
inline double trueU2_test9(int& model, double x, double y, double z) { return ay;  } // the LDG scheme in 2D
inline double trueU3_test9(int& model, double x, double y, double z) { return 0.0; } // with constant velocity field

double trueU1_test10(int& model, double x, double y, double z); // provides analytical velocity (SD test '-10')
double trueU2_test10(int& model, double x, double y, double z); // to test for convergence the LDG scheme
double trueU3_test10(int& model, double x, double y, double z); 

double trueU1_test11(int& model, double x, double y, double z); // provides analytical velocity (SD test '-11')
double trueU2_test11(int& model, double x, double y, double z); // to test for convergence the LDG scheme
double trueU3_test11(int& model, double x, double y, double z);

double trueU1_test12(int& model, double x, double y, double z); // provides analytical velocity (SD test '-12')
double trueU2_test12(int& model, double x, double y, double z); // to test for convergence the LDG scheme
double trueU3_test12(int& model, double x, double y, double z); 

double trueU1_test13(int& model, double x, double y, double z); // provides analytical velocity 
double trueU2_test13(int& model, double x, double y, double z); // (SD tests '-10', '-11', '-12')
double trueU3_test13(int& model, double x, double y, double z); // to test for convergence the LDG scheme

double trueU1_test14(int& model, double x, double y, double z); // provides analytical velocity (SD test '-100') 
double trueU2_test14(int& model, double x, double y, double z);
double trueU3_test14(int& model, double x, double y, double z);

//------------------------------------------------------------------------------------------------



//-------------------------------- assign velocity values ----------------------------------------
// these are called from transport_()
// ( used primarily for debugging )
void set_velocity(eqn *Ptr, double vel[]);
void set_velocity(eqn *Ptr, double v1, double v2, double v3);


// this is to be called from parceltst ( the main program )
//
// used in the main program to provide a velocity field for the LDG scheme
// if the flow calculations are "switched off" with velflag=0 
extern "C" {
  void setvel_( int& model, int &testnumber, int &n1, int &n2, int &n3,
                double *x1, double *x2, double *x3,
		double *u1, double *u2, double *u3 );
};
//------------------------------------------------------------------------------------------------


//------------------------------ apply the initial condition -------------------------------------
void apply_IC(int& model, int nx, double x[], int ny, double y[], int nz, double z[], double c[][8], int testnumber);
void apply_IC_discts(int& model, int nx, double x[], int ny, double y[], int nz, double z[],
                     double c[][8], int testnumber);
//------------------------------------------------------------------------------------------------


double c_circle(double& x, double& y);


#endif
