
#include "problem.h"

double c_inflow(int& model, double x, double y, double z, double t, int testnumber)
{
   double c_val;
   switch (testnumber) {
      case 1 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 2 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 3 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 4 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 5 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 6 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 7 : c_val = ( fabs(x-0.0) < 1.0e-16 ) ? 1.0 : 0.0; break;
      case 8 : c_val = 0.0; break;
      case 10: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 11: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 12: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 13: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 14: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 15: c_val = ( (y>0.65) && (y<0.75) ) ? 1.0 : 0.0; break;
      case 16: c_val = 0.0; break;
      default : c_val = 0.0;
   }
   return c_val;
}

//------------------------------------------------------------------------------------------------

double c_initial(int& model, double x, double y, double z, int testnumber)
{
   double c_val, r, r2, xx, yy, t=0.0;
   switch (testnumber) {
      case 1 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 2 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 3 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 4 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 5 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 6 : c_val = c_true(model,x,y,z,t,testnumber); break;
      case 7 : c_val = 0.0; break;
      case 8 : c_val = ((0.375<x)&&(x<0.5)) ? 1.0 : 0.0; break;
      case 10: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 11: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 12: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 13: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 14: c_val = c_true(model,x,y,z,t,testnumber); break;
      case 15: c_val = 0.0; break;
      case 16: c_val = c_circle(x,y); break;
      default: c_val = 0.0;
   }
   return c_val;   
}

//------------------------------------------------------------------------------------------------

double f(int& model, double x, double y, double z, double t, int testnumber, eqn *Ptr)
{
   double f_val, r, r2, u, DiffusionConst, D1, D2, ksi, chi, theta, omega;
   if ( (testnumber==1)||(testnumber==3)||(testnumber==5) ) { f_val = 0.0; } // purely hyperbolic tests
   else if (testnumber==2) {
      DiffusionConst = Ptr->D[0][0];
      u = trueU1_test2(model,x,y,z);
      f_val = DiffusionConst * 4.0*PI*PI * sin(2.0*PI*(x-u*t));
   }
   else if (testnumber==4) {
      DiffusionConst = Ptr->D[0][0];
      u = trueU1_test4(model,x, y, z);
      r = x - x0 - u*t; r2 = r*r;
      f_val = DiffusionConst * 2.0*K1*K2 * ( 1.0 - 2.0*K2*r2) * exp(-K2*r2); 
   }
   else if (testnumber==6) {
      DiffusionConst = Ptr->D[0][0];
      f_val = - 2.0 * DiffusionConst;
   }
   else if (testnumber==9) {
      D1 = Ptr->D[0][0]; D2 = Ptr->D[1][1];
      double fterm1 = (cos(PI*x) + cos(PI*y)) / PI;
      double fterm2 = -t*sin(PI*x)*trueU1_test9(model,x,y,z) - t*sin(PI*y)*trueU2_test9(model,x,y,z);
      double fterm3 = PI*t  *  ( D1*cos(PI*x) + D2*cos(PI*y) );
      f_val = fterm1 + fterm2 + fterm3;
   }
   else if (testnumber==10) {
      double fterm1, fterm2, fterm3;
      D1 = Ptr->D[0][0]; D2 = Ptr->D[1][1];
      fterm1 = (cos(PI*x) + cos(PI*y)) / PI;
      if (model==2) {
	 fterm2 = -t*sin(PI*x)*trueU1_test10(model,x,y,z) - t*sin(PI*y)*trueU2_test10(model,x,y,z);
      }
      else {
	 int pr=2; theta = sdparameters_(&pr); pr=3; chi = sdparameters_(&pr);
	 fterm2 = (1.0-theta)*chi * c_true(model, x, y, z, t, testnumber)
	          -t*sin(PI*x)*trueU1_test10(model,x,y,z) - t*sin(PI*y)*trueU2_test10(model,x,y,z);
      }
      fterm3 = PI*t  *  ( D1*cos(PI*x) + D2*cos(PI*y) );
      f_val = fterm1 + fterm2 + fterm3;
   }
   else if (testnumber==11) {
      double fterm1, fterm2, fterm3;
      double xx = x - t*ax, yy = y - t*ay;
      D1 = Ptr->D[0][0]; D2 = Ptr->D[1][1];
      fterm1 = (cos(PI*x) + cos(PI*y)) / PI;
      if (model==2) {
	 fterm2 = -t*sin(PI*x)*trueU1_test11(model,x,y,z) - t*sin(PI*y)*trueU2_test11(model,x,y,z);
      }
      else {
	 int pr=1; ksi = sdparameters_(&pr); pr=3; chi = sdparameters_(&pr);
	 fterm2 = (ksi - 0.5 + chi) * c_true(model, x, y, z, t, testnumber)
	          -t*sin(PI*x)*trueU1_test11(model,x,y,z) - t*sin(PI*y)*trueU2_test11(model,x,y,z);
      }
      fterm3 = PI*t  *  ( D1*cos(PI*x) + D2*cos(PI*y) );
      f_val = fterm1 + fterm2 + fterm3;
   }
   else if (testnumber==12) {
      double fterm1, fterm2, fterm3;
      D1 = Ptr->D[0][0]; D2 = Ptr->D[1][1];
      fterm1 = (cos(PI*x) + cos(PI*y)) / PI;    
      fterm2 = -t*sin(PI*x)*trueU1_test12(model,x,y,z) - t*sin(PI*y)*trueU2_test12(model,x,y,z);
      fterm3 = PI*t  *  ( D1*cos(PI*x) + D2*cos(PI*y) );
      f_val = fterm1 + fterm2 + fterm3;
   }
   else if (testnumber==14) {
      double fterm1, fterm2, fterm3;
      D1 = Ptr->D[0][0]; D2 = Ptr->D[1][1];
      fterm1 = (cos(PI*x) + cos(PI*y)) / PI;
      if (model==2) {
	 fterm2 = -t*sin(PI*x)*trueU1_test14(model,x,y,z) - t*sin(PI*y)*trueU2_test14(model,x,y,z);
      }
      else { 
	 int pr=3; chi = sdparameters_(&pr); pr=5; omega = sdparameters_(&pr);
	 fterm2 = (chi - omega*omega*sin(omega*x)*y) * c_true(model, x, y, z, t, testnumber)
	          -t*sin(PI*x)*trueU1_test14(model,x,y,z) - t*sin(PI*y)*trueU2_test14(model,x,y,z);
      }
      fterm3 = PI*t  *  ( D1*cos(PI*x) + D2*cos(PI*y) );
      f_val = fterm1 + fterm2 + fterm3;
   }
   else if ( (testnumber==13)||(testnumber==15)||(testnumber==16) ) { f_val = 0.0; }
   else { f_val = 0.0; }
   return f_val;
}

//------------------------------------------------------------------------------------------------

double c_true(int& model, double x, double y, double z, double t, int testnumber)
{
   double c_val, r, u;
   if (testnumber==1) {
      u = trueU1_test1(model,x,y,z);
      c_val = sin(2.0*PI*(x-u*t)) + 2.0;
   }
   else if (testnumber==2) {
      u = trueU1_test2(model,x,y,z);
      c_val = sin(2.0*PI*(x-u*t)) + 2.0;
   }
   else if (testnumber==3) {
      u = trueU1_test3(model,x, y, z);
      r = x - x0 - u*t;
      c_val = K1 * exp( -K2*r*r );
   }
   else if (testnumber==4) {
      u = trueU1_test4(model,x, y, z);
      r = x - x0 - u*t;
      c_val = K1 * exp( -K2*r*r );
   }
   else if (testnumber==5) {
      u = trueU1_test5(model,x, y, z);
      r = x - u * t;
      c_val = r * r;
   }
   else if (testnumber==6) {
      u = trueU1_test6(model,x, y, z);
      r = x - u * t;
      c_val = r * r;
   }
   else if ((testnumber>=9)&&(testnumber<=14)) {
      c_val =  t *  ( cos(PI*x) + cos(PI*y) )  / PI;
   }
   else { c_val = 0.0; }
   return c_val;
}

double qx_true(int& model, eqn *Ptr, double x, double y, double z, double t, int testnumber)
{
   double qx, r, u, DiffusionConst;
   if ( (testnumber==1)||(testnumber==3)||(testnumber==5) ) {
      qx = 0.0;
   }
   else if (testnumber==2) {
      DiffusionConst = Ptr->D[0][0];
      u = trueU1_test2(model,x,y,z);
      qx = -DiffusionConst * 2.0*PI * cos(2.0*PI*(x-u*t)); 
   }
   else if (testnumber==4) {
      DiffusionConst = Ptr->D[0][0];
      u = trueU1_test4(model,x, y, z);
      r = x - x0 - u*t;
      qx = DiffusionConst * 2.0*K1*K2*r*exp(-K2*r*r);
   }
   else if (testnumber==6) {
      DiffusionConst = Ptr->D[0][0];
      u = trueU1_test6(model,x, y, z);
      qx = -DiffusionConst * 2.0 * (x - u * t); 
   }
   else if (((testnumber>=9)&&(testnumber<=12))||(testnumber==14)) {
      DiffusionConst = Ptr->D[0][0];
      qx = DiffusionConst * sin(PI*x) * t;
   }
   else { qx = 0.0; }
   return qx;
}

double qy_true(int& model, eqn *Ptr, double x, double y, double z, double t, int testnumber)
{
   double qy, DiffusionConst;
   if ( (testnumber>=1)&&(testnumber<=6) ) {
      qy = 0.0;
   }
   else if (((testnumber>=9)&&(testnumber<=12))||(testnumber==14)) {
      DiffusionConst = Ptr->D[1][1];
      qy = DiffusionConst * sin(PI*y) * t;
   }
   else { qy = 0.0; }
   return qy;
}

double qz_true(int& model, eqn *Ptr, double x, double y, double z, double t, int testnumber)
{
   double qz, DiffusionConst;
   if ( (testnumber>=1)&&(testnumber<=6) ) {
      qz = 0.0;
   }
   else if (((testnumber>=9)&&(testnumber<=12))||(testnumber==14)) {
      qz = 0.0;
   }
   else { qz = 0.0; }
   return qz;
}

//------------------------------------------------------------------------------------------------

double trueU1_test10(int& model, double x, double y, double z)
{
   double u1, ksi, beta, chi;
   int pr;
   if (model==2) {
      pr=1; ksi = sdparameters_(&pr);
      u1 = (2.0-x)*(1.5-y)*(y-ksi);
   }
   else {
      pr = 2; beta = sdparameters_(&pr);
      pr = 3; chi  = sdparameters_(&pr);
      u1 = beta*chi*(2.0-x); 
   }
   return u1;
}
double trueU2_test10(int& model, double x, double y, double z)
{
   double u2, ksi, chi;
   int pr;
   if (model==2) {
      pr=1; ksi = sdparameters_(&pr);
      double yy = y*y;
      u2 = - y*yy/3.0 + (ksi+1.5)*yy/2.0 - 1.5*ksi*y - 0.5;
   }
   else {
      pr = 3; chi  = sdparameters_(&pr);
      u2 = chi*(y+0.5); 
   }
   return u2;
}
double trueU3_test10(int& model, double x, double y, double z) { return 0.0; }
//------------------------------------------------------------------------------------------------
double trueU1_test11(int& model, double x, double y, double z)
{
   double u1, ksi, beta, chi;
   int pr;
   if (model==2) {
      pr=1; ksi = sdparameters_(&pr);
      u1 = (2.0-x)*(1.5-y)*(y-ksi);
   }
   else {
      pr=1; ksi = sdparameters_(&pr);
      u1 = (2.0-x)*(0.5-ksi); 
   }
   return u1;
}
double trueU2_test11(int& model, double x, double y, double z)
{
   double u2, ksi, chi;
   int pr;
   if (model==2) {
      pr=1; ksi = sdparameters_(&pr);
      double yy = y*y;
      u2 = - y*yy/3.0 + (ksi+1.5)*yy/2.0 - 1.5*ksi*y - 0.5;
   }
   else {
      pr = 3; chi  = sdparameters_(&pr);
      u2 = chi*(y+0.5); 
   }
   return u2;
}
double trueU3_test11(int& model, double x, double y, double z) { return 0.0; }

//------------------------------------------------------------------------------------------------
double trueU1_test12(int& model, double x, double y, double z)
{
   double u1, omega, BJScoeff;
   int pr;
   
   pr=4; BJScoeff = sdparameters_(&pr);
   pr=5; omega    = sdparameters_(&pr);

   u1 = sin(x/BJScoeff + omega)*exp(y/BJScoeff);
   return u1;
}
double trueU2_test12(int& model, double x, double y, double z)
{
   double u2, omega, BJScoeff;
   int pr;
   
   pr=4; BJScoeff = sdparameters_(&pr);
   pr=5; omega    = sdparameters_(&pr);

   u2 = -cos(x/BJScoeff + omega)*exp(y/BJScoeff);      
   return u2;
}
double trueU3_test12(int& model, double x, double y, double z) { return 0.0; }

//------------------------------------------------------------------------------------------------

double trueU1_test13(int& model, double x, double y, double z) { return truesdvel1_(&model,&x,&y,&z); }
double trueU2_test13(int& model, double x, double y, double z) { return truesdvel2_(&model,&x,&y,&z); }
double trueU3_test13(int& model, double x, double y, double z) { return 0.0; }

//------------------------------------------------------------------------------------------------

double trueU1_test14(int& model, double x, double y, double z)
{
   double u1, ksi, omega;
   int pr;
   if (model==2) {
      pr=1; ksi = sdparameters_(&pr);
      u1 = (2.0-x)*(1.5-y)*(y-ksi);
   }
   else {
      pr=5; omega = sdparameters_(&pr);
      u1 = omega*cos(omega*x)*y;
   }
   return u1;
}

double trueU2_test14(int& model, double x, double y, double z)
{
   double u2, ksi, chi, omega;
   int pr;
   if (model==2) {
      double y2 = y*y, y3 = y2*y;
      pr=1; ksi = sdparameters_(&pr);
      pr=5; omega = sdparameters_(&pr);
      u2 = (-y3)/3.0 + (1.5+ksi)*(y2)/2.0 - 1.5*ksi*y - 0.5 + sin(omega*x);
   }
   else {
      pr=3; chi   = sdparameters_(&pr);
      pr=5; omega = sdparameters_(&pr);
      u2 = chi * (y + 0.5) + sin(omega*x);
   }
   return u2;
}

double trueU3_test14(int& model, double x, double y, double z) { return 0.0; }

//------------------------------------------------------------------------------------------------

void set_velocity(eqn *Ptr, double v1, double v2, double v3)
{
   int i, n=(Ptr->n1) * (Ptr->n2) * (Ptr->n3);
   for (i=0; i<(8*n); ++i) Ptr->u1[i] = v1;
   for (i=0; i<(8*n); ++i) Ptr->u2[i] = v2;
   for (i=0; i<(8*n); ++i) Ptr->u3[i] = v3;
}

void set_velocity(eqn *Ptr, double vel[])
{
   int i, n=(Ptr->n1) * (Ptr->n2) * (Ptr->n3);
   for (i=0; i<(8*n); ++i) Ptr->u1[i] = vel[0];
   for (i=0; i<(8*n); ++i) Ptr->u2[i] = vel[1];
   for (i=0; i<(8*n); ++i) Ptr->u3[i] = vel[2];
}

void setvel_( int& model, int &testnumber, int &n1, int &n2, int &n3,
                    double *x1, double *x2, double *x3,
		    double *u1, double *u2, double *u3 )
{

   velfun Vel1[velN], Vel2[velN], Vel3[velN];
   Vel1[0]  = trueU1_test1;  Vel2[0]  = trueU2_test1;  Vel3[0]  = trueU3_test1;
   Vel1[1]  = trueU1_test2;  Vel2[1]  = trueU2_test2;  Vel3[1]  = trueU3_test2;
   Vel1[2]  = trueU1_test3;  Vel2[2]  = trueU2_test3;  Vel3[2]  = trueU3_test3;
   Vel1[3]  = trueU1_test4;  Vel2[3]  = trueU2_test4;  Vel3[3]  = trueU3_test4;
   Vel1[4]  = trueU1_test5;  Vel2[4]  = trueU2_test5;  Vel3[4]  = trueU3_test5;
   Vel1[5]  = trueU1_test6;  Vel2[5]  = trueU2_test6;  Vel3[5]  = trueU3_test6;
   Vel1[6]  = trueU1_test7;  Vel2[6]  = trueU2_test7;  Vel3[6]  = trueU3_test7;
   Vel1[7]  = trueU1_test8;  Vel2[7]  = trueU2_test8;  Vel3[7]  = trueU3_test8;
   Vel1[8]  = trueU1_test9;  Vel2[8]  = trueU2_test9;  Vel3[8]  = trueU3_test9;
   Vel1[9]  = trueU1_test10; Vel2[9]  = trueU2_test10; Vel3[9]  = trueU3_test10;
   Vel1[10] = trueU1_test11; Vel2[10] = trueU2_test11; Vel3[10] = trueU3_test11;
   Vel1[11] = trueU1_test12; Vel2[11] = trueU2_test12; Vel3[11] = trueU3_test12;
   Vel1[12] = trueU1_test13; Vel2[12] = trueU2_test13; Vel3[12] = trueU3_test13;
   Vel1[13] = trueU1_test14; Vel2[13] = trueU2_test14; Vel3[13] = trueU3_test14;

   if( (testnumber>velN)||(testnumber<1) ) {
      cout << "Transport test case " << testnumber << " does not have analytical velocity field available!" << endl;
      exit(1);
   }

   int t = testnumber - 1; // indexing of the arrays starts from 0
   int n = n1*n2*n3;
   for (int i3=0; i3<n3; ++i3 ) {
      for (int i2=0; i2<n2; ++i2 ) {
         for (int i1=0; i1<n1; ++i1 ) {

            int elindex = i1 + n1*i2 + n1*n2*i3;

            u1[elindex+0*n] = Vel1[t]( model, x1[i1],   x2[i2],   x3[i3]   );
            u1[elindex+1*n] = Vel1[t]( model, x1[i1+1], x2[i2],   x3[i3]   );
            u1[elindex+2*n] = Vel1[t]( model, x1[i1+1], x2[i2+1], x3[i3]   );
            u1[elindex+3*n] = Vel1[t]( model, x1[i1],   x2[i2+1], x3[i3]   );
            u1[elindex+4*n] = Vel1[t]( model, x1[i1],   x2[i2],   x3[i3+1] );
            u1[elindex+5*n] = Vel1[t]( model, x1[i1+1], x2[i2],   x3[i3+1] );
            u1[elindex+6*n] = Vel1[t]( model, x1[i1+1], x2[i2+1], x3[i3+1] );
            u1[elindex+7*n] = Vel1[t]( model, x1[i1],   x2[i2+1], x3[i3+1] );

	    u2[elindex+0*n] = Vel2[t]( model, x1[i1],   x2[i2],   x3[i3]   );
            u2[elindex+1*n] = Vel2[t]( model, x1[i1+1], x2[i2],   x3[i3]   );
            u2[elindex+2*n] = Vel2[t]( model, x1[i1+1], x2[i2+1], x3[i3]   );
            u2[elindex+3*n] = Vel2[t]( model, x1[i1],   x2[i2+1], x3[i3]   );
            u2[elindex+4*n] = Vel2[t]( model, x1[i1],   x2[i2],   x3[i3+1] );
            u2[elindex+5*n] = Vel2[t]( model, x1[i1+1], x2[i2],   x3[i3+1] );
            u2[elindex+6*n] = Vel2[t]( model, x1[i1+1], x2[i2+1], x3[i3+1] );
            u2[elindex+7*n] = Vel2[t]( model, x1[i1],   x2[i2+1], x3[i3+1] );

	    u3[elindex+0*n] = Vel3[t]( model, x1[i1],   x2[i2],   x3[i3]   );
            u3[elindex+1*n] = Vel3[t]( model, x1[i1+1], x2[i2],   x3[i3]   );
            u3[elindex+2*n] = Vel3[t]( model, x1[i1+1], x2[i2+1], x3[i3]   );
            u3[elindex+3*n] = Vel3[t]( model, x1[i1],   x2[i2+1], x3[i3]   );
            u3[elindex+4*n] = Vel3[t]( model, x1[i1],   x2[i2],   x3[i3+1] );
            u3[elindex+5*n] = Vel3[t]( model, x1[i1+1], x2[i2],   x3[i3+1] );
            u3[elindex+6*n] = Vel3[t]( model, x1[i1+1], x2[i2+1], x3[i3+1] );
            u3[elindex+7*n] = Vel3[t]( model, x1[i1],   x2[i2+1], x3[i3+1] );

	 }
      }
   }
	       
}

//------------------------------------------------------------------------------------------------

void apply_IC(int& model, int nx, double x[], int ny, double y[], int nz, double z[], double c[][8], int testnumber)
// n?     :   number of divisions in the local subdomain along the corresponding coordinate axis
// x,y,z  :   arrays of length (n? + 1) to hold the coordinates of the nodes
{
   for (int iz=0; iz<nz; ++iz ) {
      for (int iy=0; iy<ny; ++iy ) {
         for (int ix=0; ix<nx; ++ix ) {
            int elindex = ix + nx*iy + nx*ny*iz;

            // the local ordering is explained in utilities.h
            
            //cout << "element " << elindex << endl;
            c[elindex][0] = c_initial( model, x[ix],   y[iy],   z[iz],   testnumber ); //cout << c[elindex][0] << " ";
            c[elindex][1] = c_initial( model, x[ix+1], y[iy],   z[iz],   testnumber ); //cout << c[elindex][1] << " ";
            c[elindex][2] = c_initial( model, x[ix+1], y[iy+1], z[iz],   testnumber ); //cout << c[elindex][2] << " ";
            c[elindex][3] = c_initial( model, x[ix],   y[iy+1], z[iz],   testnumber ); //cout << c[elindex][3] << " ";
            c[elindex][4] = c_initial( model, x[ix],   y[iy],   z[iz+1], testnumber ); //cout << c[elindex][4] << " ";
            c[elindex][5] = c_initial( model, x[ix+1], y[iy],   z[iz+1], testnumber ); //cout << c[elindex][5] << " ";
            c[elindex][6] = c_initial( model, x[ix+1], y[iy+1], z[iz+1], testnumber ); //cout << c[elindex][6] << " ";
            c[elindex][7] = c_initial( model, x[ix],   y[iy+1], z[iz+1], testnumber ); //cout << c[elindex][7] << " " << endl;
            //cout << "***********" << endl;
         }
      }
   }
}

//------------------------------------------------------------------------------------------------

void apply_IC_discts(int& model, int nx, double x[], int ny, double y[], int nz, double z[],
                     double c[][8], int testnumber)
{
   const double EPS = 1.0e-10;

   for (int iz=0; iz<nz; ++iz ) {
      for (int iy=0; iy<ny; ++iy ) {
         for (int ix=0; ix<nx; ++ix ) {
            int elindex = ix + nx*iy + nx*ny*iz;

            c[elindex][0] = c_initial( model, x[ix]  +EPS,  y[iy]  +EPS,  z[iz]  +EPS, testnumber );
            c[elindex][1] = c_initial( model, x[ix+1]-EPS,  y[iy]  +EPS,  z[iz]  +EPS, testnumber );
            c[elindex][2] = c_initial( model, x[ix+1]-EPS,  y[iy+1]-EPS,  z[iz]  +EPS, testnumber );
            c[elindex][3] = c_initial( model, x[ix]  +EPS,  y[iy+1]-EPS,  z[iz]  +EPS, testnumber );
            c[elindex][4] = c_initial( model, x[ix]  +EPS,  y[iy]  +EPS,  z[iz+1]-EPS, testnumber );
            c[elindex][5] = c_initial( model, x[ix+1]-EPS,  y[iy]  +EPS,  z[iz+1]-EPS, testnumber );
            c[elindex][6] = c_initial( model, x[ix+1]-EPS,  y[iy+1]-EPS,  z[iz+1]-EPS, testnumber );
            c[elindex][7] = c_initial( model, x[ix]  +EPS,  y[iy+1]-EPS,  z[iz+1]-EPS, testnumber );

         }
      }
   }   
}

//------------------------------------------------------------------------------------------------

double c_circle(double& x, double& y)
{
   double X = (x-0.2);
   double Y = (y-0.7);
   double R = sqrt(X*X + Y*Y);
   double c_val = (R<0.1) ? 1.0 : 0.0;
   return c_val;
}
