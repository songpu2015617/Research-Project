#ifndef _QUAD_H
#define _QUAD_H

#include "common.h"
#include "basis.h"


//------------------------------------------------------------------------------------------------
//                           ^
// compute (wi, wj)    where E is the reference element
//                 ^
//                 E 
double quad(function wi, function wj);
//------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------
//                           ^
// compute (wi, wj)    where F  =  [-1, 1] x [-1, 1]
//                 ^
//                 F
// choose = 0, 1 or 2 and represents the coordinate being fixed to be equal to val (x, y or z respectively) 
double quad(int choose, double val, function wi, function wj);
//------------------------------------------------------------------------------------------------


#endif
