#ifndef _LDG_H
#define _LDG_H

#include "utilities.h"
#include "limiter.h"


//------------------------------------------------------------------------------------------------
// compute (wi ,wj )   where E is an element
//                  E
double quad(int& i1, int& i2, int& i3, function wi, function wj, eqn *Ptr);
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
// compute   term1  =   (  C U + Q , grad w  )    for i = 0, 1, ..., 7
//                                         i  E

void LDG_term1(int& i1, int& i2, int& i3, eqn *Ptr,
               double c[][8], double qx[][8], double qy[][8], double qz[][8], double term[8]);

//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
//                               ^    _
// compute   term2  =   <  C U . n , w  >    __        for i = 0, 1, ..., 7
//                                    i  dE |  | Gamma
//                            -
// Note: C here is is either C or C depending on the boundary type
//                                 I
//
void LDG_term2(int& model, int& i1, int& i2, int& i3, eqn *Ptr, int *mynbr, 
               double c[][8], double term[8], int& testnumber, double& CurrentTime);

//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
//                            upwd       _      ^    _
// compute   term3  =   <  ( C     U  +  Q  ) . n , w   >         for i = 0, 1, ..., 7
//                                                   i   dE\Gamma
//
void LDG_term3(int& i1, int& i2, int& i3, eqn *Ptr, int *mynbr,
               double qbarx[][8], double qbary[][8], double qbarz[][8],
               double upwd1[][8], double upwd2[][8], double upwd3[][8],
               double c[][8], double term[8]);

//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
// compute   term4  =   (  f , w  )    for i = 0, 1, ..., 7
//                              i  E
//
void LDG_term4(int& model, int& i1, int& i2, int& i3,
               eqn *Ptr, double term[8], int& testnumber, double& CurrentTime);

//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
// NOTE: The integrals are computed (exactly) using the Simpson's rule
void LDG_rhs(int& model, int& i1, int& i2, int& i3, eqn *Ptr, int *mynbr, double c[][8],
             double upwd1[][8], double upwd2[][8], double upwd3[][8],
             double qx[][8], double qy[][8], double qz[][8],
	     double qbarx[][8], double qbary[][8], double qbarz[][8],
             double TimeStep, double rhs[8], int& testnumber, double& CurrentTime);
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
void fwd_euler( int& mynod, int& model, int *mynbr, eqn *Ptr, int& limflag, double& alpha, int& dflag, int& dtype,
                double c[][8], double cexp[][8], double cvar[][8], double TimeStep, double TimeMax, int skip,
                fstream& out, int& testnumber, int& ldgerr, int& StochFlag, int& StochCount, double& cweight,
                int *RealizationIdx, int *NumOfRealizations, int *totalloops );
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
void rk2( int& mynod, int& model, int *mynbr, eqn *Ptr, int& limflag, double& alpha, int& dflag, int& dtype,
          double c[][8], double cexp[][8], double cvar[][8], double TimeStep, double TimeMax, int skip,
	      fstream& out, int& testnumber, int& ldgerr, int& StochFlag, int& StochCount, double& cweight,
          int *RealizationIdx, int *NumOfRealizations, int *totalloops );
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
void rk4( int& mynod, int& model, int *mynbr, eqn *Ptr, int& limflag, double& alpha, int& dflag, int& dtype,
          double c[][8], double cexp[][8], double cvar[][8], double TimeStep, double TimeMax, int skip,
	      fstream& out, int& testnumber, int& ldgerr, int& StochFlag, int& StochCount, double& cweight,
          int *RealizationIdx, int *NumOfRealizations, int *totalloops );
//------------------------------------------------------------------------------------------------


#endif
