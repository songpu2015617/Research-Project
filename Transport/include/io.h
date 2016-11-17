#ifndef _IO_H
#define _IO_H

#include "problem.h"


//----------------------------------- read the input file ----------------------------------------
void read_file(ifstream& in, double bdry[3][2],
               int& tflag, int& velflag, int& testnumber,
               int& limflag, double& alpha, int& ldgerr, int& dflag, int& dtype, double D[],
               double& d_molecular, double& d_longitudinal, double& d_transverse,
               double& porosity, double& d_visc,
               int& IntegrationMethod, double& T, double& dt, int& skip);
//------------------------------------------------------------------------------------------------


//----------------- read the input file from parceltst ( the main program ) ----------------------
extern "C" {
  void readtparam_(int& tflag, int& velflag, int& testnumber, int& dtype);
};
//------------------------------------------------------------------------------------------------

//-------------------------- read debug data from the input file ---------------------------------
void read_file_dbg(ifstream& in, double vel[]);
//------------------------------------------------------------------------------------------------


//---------------------- print the output to a file in Tecplot format ----------------------------
void frame_name(int mynod, int current_frame_number, char frame[]); // make a name for the current frame
void realization_frame_name(int mynod, int StochCount, int current_frame_number, char frame[]);
void print_output( int& model, eqn *Ptr, double c[][8], fstream& out, int& frame, double& t, int& testnumber, int& mynod);
void print_output( eqn *Ptr, double qx[][8], double qy[][8], double qz[][8], fstream& out );

// miz17, bag8 : io functions for expectation and variance
void write_binary( eqn *Ptr, double cexp[][8], double cvar[][8], fstream& out);
void read_add_stochastic( eqn *Ptr, double c[][8], double cexp[][8], double cvar[][8],
                          double &cweight, fstream& in, bool first_realization, int &mynod);
void subtract_square_exp( eqn *Ptr, double cexp[][8], double cvar[][8]);
void print_stochastic( eqn *Ptr, double cexp[][8], double cvar[][8], fstream& out, int& frame, double& t, int& mynod );
//------------------------------------------------------------------------------------------------

#endif

