#include "quad.h"


double quad(function wi, function wj)
{
   double crd[3] = {-1.0, 0.0, 1.0};
   double weight[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0};
   double integral = 0.0;

   // loop over the quadrature points in the element
   for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
	 for(int k=0; k<3; ++k)
	    integral += weight[i]  *  weight[j]  *  weight[k]  * 
	                wi(crd[i], crd[j], crd[k])  *  wj(crd[i], crd[j], crd[k]);

   return integral;

}


//------------------------------------------------------------------------------------------------


double quad(int choose, double val, function wi, function wj)
{
   double crd[3] = {-1.0, 0.0, 1.0};
   double weight[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0};
   double integral = 0.0;
      
   switch (choose) {

   case 0:
      for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
	 integral += weight[i]  *  weight[j]  *  wi(val, crd[i], crd[j])  *  wj(val, crd[i], crd[j]);	     
      break;

   case 1:
      for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
	 integral += weight[i]  *  weight[j]  *  wi(crd[i], val, crd[j])  *  wj(crd[i], val, crd[j]);
      break;

   case 2:
      for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
	 integral += weight[i]  *  weight[j]  *  wi(crd[i], crd[j], val)  *  wj(crd[i], crd[j], val);
      break;

   }

   return integral;
}
