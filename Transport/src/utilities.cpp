#include "utilities.h"


void cholesky_solve(int n, double U[], double x[], double b[])
{
   double y[n];

   y[0] = b[0] / U[0];
   for (int i=1; i<n; ++i) {
      y[i] = b[i];
      for (int j=0; j<i; ++j)
	 y[i] -= U[j+i*n] * y[j];
      y[i] = y[i] / U[i+i*n];
   }

   x[n-1] = y[n-1] / U[(n-1) + (n-1)*n];
   for (int i=(n-2); i>=0; --i) {
      x[i] = y[i];
      for (int j=i+1; j<n; ++j)
	 x[i] -= U[i+j*n] * x[j];
      x[i] = x[i] / U[i+i*n];
   }
}



void cramer(double A[][3], double b[], double x[])
{
   int i;
   double T[3][3];
   double determinant = det(A);

   if (fabs(determinant) < 1.0E-15) {
      cout << "Encountered zero determinant in" << "\nvoid cramer(double[][3], b[], x[])" << endl;
      exit(1);
   }
   else {
      for (i=0; i<3; ++i) T[i][0] = b[i];
      for (i=0; i<3; ++i) T[i][1] = A[i][1];
      for (i=0; i<3; ++i) T[i][2] = A[i][2];
      x[0] = det(T) / determinant;

      for (i=0; i<3; ++i) T[i][0] = A[i][0];
      for (i=0; i<3; ++i) T[i][1] = b[i];
      x[1] = det(T) / determinant;

      for (i=0; i<3; ++i) T[i][1] = A[i][1];
      for (i=0; i<3; ++i) T[i][2] = b[i];
      x[2] = det(T) / determinant;
   }
}



//------------------------------------------------------------------------------------------------


void avr_x(double c[][8], double cbarx[][8], int nx, int ny, int nz)
{
   for (int iz=0; iz<nz; ++iz ) {
      for (int iy=0; iy<ny; ++iy ) {
         for (int ix=0; ix<nx; ++ix ) {

            int k = ix + nx*iy + nx*ny*iz;

            int kE = k + 1;
            int kW = k - 1;

            int W = (ix==0)      ? 0 : 1;
            int E = (ix==(nx-1)) ? 0 : 1;

            // local node 0
            cbarx[k][0] = (W==0) ? c[k][0] : (c[k][0] + c[kW][1])/2.0;
            //---------------------------------------------

            // local node 1
	    cbarx[k][1] = (E==0) ? c[k][1] : (c[k][1] + c[kE][0])/2.0;
            //---------------------------------------------

            // local node 2
	    cbarx[k][2] = (E==0) ? c[k][2] : (c[k][2] + c[kE][3])/2.0;
            //---------------------------------------------

            // local node 3
	    cbarx[k][3] = (W==0) ? c[k][3] : (c[k][3] + c[kW][2])/2.0;
            //---------------------------------------------

            // local node 4
            cbarx[k][4] = (W==0) ? c[k][4] : (c[k][4] + c[kW][5])/2.0;
            //---------------------------------------------

            // local node 5
	    cbarx[k][5] = (E==0) ? c[k][5] : (c[k][5] + c[kE][4])/2.0;
            //---------------------------------------------

            // local node 6
            cbarx[k][6] = (E==0) ? c[k][6] : (c[k][6] + c[kE][7])/2.0;
            //---------------------------------------------

            // local node 7
            cbarx[k][7] = (W==0) ? c[k][7] : (c[k][7] + c[kW][6])/2.0;
            //---------------------------------------------

         }
      }
   }
}


//------------------------------------------------------------------------------------------------


void avr_y(double c[][8], double cbary[][8], int nx, int ny, int nz)
{
   for (int iz=0; iz<nz; ++iz ) {
      for (int iy=0; iy<ny; ++iy ) {
         for (int ix=0; ix<nx; ++ix ) {

            int k = ix + nx*iy + nx*ny*iz;

            int kN = k + nx;
            int kS = k - nx;

            int S = (iy==0)      ? 0 : 1;
            int N = (iy==(ny-1)) ? 0 : 1;

            // local node 0
            cbary[k][0] = (S==0) ? c[k][0] : (c[k][0] + c[kS][3])/2.0;
            //---------------------------------------------

            // local node 1
            cbary[k][1] = (S==0) ? c[k][1] : (c[k][1] + c[kS][2])/2.0;
            //---------------------------------------------

            // local node 2
            cbary[k][2] = (N==0) ? c[k][2] : (c[k][2] + c[kN][1])/2.0;
            //---------------------------------------------

            // local node 3
	    cbary[k][3] = (N==0) ? c[k][3] : (c[k][3] + c[kN][0])/2.0;
            //---------------------------------------------

            // local node 4
            cbary[k][4] = (S==0) ? c[k][4] : (c[k][4] + c[kS][7])/2.0;
            //---------------------------------------------

            // local node 5
            cbary[k][5] = (S==0) ? c[k][5] : (c[k][5] + c[kS][6])/2.0;
            //---------------------------------------------

            // local node 6
            cbary[k][6] = (N==0) ? c[k][6] : (c[k][6] + c[kN][5])/2.0;
            //---------------------------------------------

            // local node 7
            cbary[k][7] = (N==0) ? c[k][7] : (c[k][7] + c[kN][4])/2.0;
            //---------------------------------------------

         }
      }
   }
}


//------------------------------------------------------------------------------------------------


void avr_z(double c[][8], double cbarz[][8], int nx, int ny, int nz)
{
   for (int iz=0; iz<nz; ++iz ) {
      for (int iy=0; iy<ny; ++iy ) {
         for (int ix=0; ix<nx; ++ix ) {

            int k = ix + nx*iy + nx*ny*iz;

            int kU = k + nx*ny;
            int kD = k - nx*ny;

            int D = (iz==0)      ? 0 : 1;
            int U = (iz==(nz-1)) ? 0 : 1;

            // local node 0
            cbarz[k][0] = (D==0) ? c[k][0] : (c[k][0] + c[kD][4])/2.0;
            //---------------------------------------------

            // local node 1
            cbarz[k][1] = (D==0) ? c[k][1] : (c[k][1] + c[kD][5])/2.0;
            //---------------------------------------------

            // local node 2
            cbarz[k][2] = (D==0) ? c[k][2] : (c[k][2] + c[kD][6])/2.0;
            //---------------------------------------------

            // local node 3
            cbarz[k][3] = (D==0) ? c[k][3] : (c[k][3] + c[kD][7])/2.0;
            //---------------------------------------------

            // local node 4
            cbarz[k][4] = (U==0) ? c[k][4] : (c[k][4] + c[kU][0])/2.0;
            //---------------------------------------------

            // local node 5
            cbarz[k][5] = (U==0) ? c[k][5] : (c[k][5] + c[kU][1])/2.0;
            //---------------------------------------------

            // local node 6
            cbarz[k][6] = (U==0) ? c[k][6] : (c[k][6] + c[kU][2])/2.0;
            //---------------------------------------------

            // local node 7
            cbarz[k][7] = (U==0) ? c[k][7] : (c[k][7] + c[kU][3])/2.0;
            //---------------------------------------------

         }
      }
   }
}


//------------------------------------------------------------------------------------------------


void load_buffers(double c[][8], double bf1[], double bf2[], double bf3[], int nx, int ny, int nz)
{
   int k, ix, iy, iz, count;

   ix = 0, count = 0;
   for (iz=0; iz<nz; ++iz ) {
      for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         bf1[count] = c[k][0]; ++count; bf1[count] = c[k][3]; ++count;
      }
      for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         bf1[count] = c[k][4]; ++count; bf1[count] = c[k][7]; ++count;
      }
   }
   ix = nx-1;
   for (iz=0; iz<nz; ++iz ) {
      for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         bf1[count] = c[k][1]; ++count; bf1[count] = c[k][2]; ++count;
       }
       for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         bf1[count] = c[k][5]; ++count; bf1[count] = c[k][6]; ++count;
       }
   }

   //----------------------------------

   iy = 0; count = 0;
   for (iz=0; iz<nz; ++iz ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf2[count] = c[k][0]; ++count; bf2[count] = c[k][1]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf2[count] = c[k][4]; ++count; bf2[count] = c[k][5]; ++count;
      }
   }
   iy = ny-1;
   for (iz=0; iz<nz; ++iz ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf2[count] = c[k][3]; ++count; bf2[count] = c[k][2]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf2[count] = c[k][7]; ++count; bf2[count] = c[k][6]; ++count;
      }
   }

   //----------------------------------

   iz = 0; count =0;
   for (iy=0; iy<ny; ++iy ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf3[count] = c[k][0]; ++count; bf3[count] = c[k][1]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf3[count] = c[k][3]; ++count; bf3[count] = c[k][2]; ++count;
      }
   }
   iz = nz-1;
   for (iy=0; iy<ny; ++iy ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf3[count] = c[k][4]; ++count; bf3[count] = c[k][5]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         bf3[count] = c[k][7]; ++count; bf3[count] = c[k][6]; ++count;
      }
   }

}


//------------------------------------------------------------------------------------------------


void read_buffers(double cbarx[][8], double cbary[][8], double cbarz[][8],
                    double bf1[], double bf2[], double bf3[],
                    int nx, int ny, int nz)
{
   int k, ix, iy, iz, count;

   ix = 0; count = 0;
   for (iz=0; iz<nz; ++iz ) {
      for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarx[k][0] = bf1[count]; ++count; cbarx[k][3] = bf1[count]; ++count;
      }
      for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarx[k][4] = bf1[count]; ++count; cbarx[k][7] = bf1[count]; ++count;
      }
   }
   ix = nx-1;
   for (iz=0; iz<nz; ++iz ) {
      for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarx[k][1] = bf1[count]; ++count; cbarx[k][2] = bf1[count]; ++count;
      }
      for (iy=0; iy<ny; ++iy ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarx[k][5] = bf1[count]; ++count; cbarx[k][6] = bf1[count]; ++count;
      }
   }

   //----------------------------------

   iy = 0; count = 0;
   for (iz=0; iz<nz; ++iz ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbary[k][0] = bf2[count]; ++count; cbary[k][1] = bf2[count]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbary[k][4] = bf2[count]; ++count; cbary[k][5] = bf2[count]; ++count;
      }
   }
   iy = ny-1;
   for (iz=0; iz<nz; ++iz ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbary[k][3] = bf2[count]; ++count; cbary[k][2] = bf2[count]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbary[k][7] = bf2[count]; ++count; cbary[k][6] = bf2[count]; ++count;
      }
   }

   //----------------------------------

   iz = 0; count =0;
   for (iy=0; iy<ny; ++iy ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarz[k][0] = bf3[count]; ++count; cbarz[k][1] = bf3[count]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarz[k][3] = bf3[count]; ++count; cbarz[k][2] = bf3[count]; ++count;
      }
   }
   iz = nz-1;
   for (iy=0; iy<ny; ++iy ) {
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarz[k][4] = bf3[count]; ++count; cbarz[k][5] = bf3[count]; ++count;
      }
      for (ix=0; ix<nx; ++ix ) {
         k = ix + nx*iy + nx*ny*iz;
         cbarz[k][7] = bf3[count]; ++count; cbarz[k][6] = bf3[count]; ++count;
      }
   }

}


//------------------------------------------------------------------------------------------------


void swapbdry( int n1, int n2, int n3, int *mynbr, int mynbrSize[],
               double *bf1, double *bf2, double *bf3, double *df1, double *df2, double *df3,
               int BdryVarType, int AvgType)
{
   if (AvgType!=1) { copyiface_(mynbr,n1,n2,n3,bf1,bf2,bf3,bf1,bf2,bf3); }

   int n12 = n1*n2, n13 = n1*n3, n23 = n2*n3, s;
   s = ( n12 > n13 ) ? n12 : n13;
   s = 16  *  (  (s > n23) ? s : n23  );

   double *a1 = new double [s];
   double *a2 = new double [s];

   double *dptr = NULL;
   int *iptr = NULL;

   // The value of FLUXTYPE here must be the same as in /Include/msgtypes.hf
   int FLUXTYPE = 43;
   int mortarflag = 0;
   int mortartype[6] = {0,0,0,0,0,0};

   swapface_(mynbr,mynbrSize,n1,n2,n3,bf1,bf2,bf3,df1,df2,df3,s,a1,a2,FLUXTYPE,BdryVarType,
             dptr,dptr,dptr,mortarflag,mortartype,iptr,dptr,dptr,dptr,dptr,dptr,dptr,
             dptr,dptr,dptr,dptr,dptr,dptr,dptr,dptr,dptr,dptr,dptr,dptr,AvgType,
             iptr,iptr,iptr,iptr,iptr,iptr,iptr,iptr,iptr,iptr,iptr,iptr,
	         iptr,iptr,iptr,iptr,iptr,iptr,dptr,dptr,dptr,dptr,dptr,dptr);

   delete a2;
   delete a1;

      
}


//------------------------------------------------------------------------------------------------


void interface_val( double c[][8], double upwd1[][8], double upwd2[][8], double upwd3[][8],
                    int n1, int n2, int n3, int *mynbr, fstream& out )
{
   // needed by swapbdry()
   double *dptr = NULL;
   int *iptr = NULL;
   int mortar=0, AvgType=4, BdryVarType=2;
   int mynbrSize[12];

   int nodesX = 2*n1, nodesY = 2*n2, nodesZ = 2*n3; // # of nodes in one direction
   getnbrsize_(mynbr, mynbrSize, nodesX, nodesY, nodesZ);

   double (*bf1) = new double [2*nodesY*nodesZ];
   double (*bf2) = new double [2*nodesX*nodesZ];
   double (*bf3) = new double [2*nodesX*nodesY];
   double (*df1) = new double [2*nodesY*nodesZ];
   double (*df2) = new double [2*nodesX*nodesZ];
   double (*df3) = new double [2*nodesX*nodesY];
   
   load_buffers(c, bf1, bf2, bf3, n1, n2, n3);
   
   swapbdry( nodesX, nodesY, nodesZ, mynbr, mynbrSize, bf1, bf2, bf3, df1, df2, df3, BdryVarType, AvgType);

   int k, i1, i2, i3, count;

   count = 0;
   for (i3=0; i3<n3; ++i3 ) {
      for (i2=0; i2<n2; ++i2 ) {
         k = i2 + n2*i3;
         upwd1[k][0] = df1[count]; ++count; upwd1[k][3] = df1[count]; ++count;
      }
      for (i2=0; i2<n2; ++i2 ) {
         k = i2 + n2*i3;
         upwd1[k][4] = df1[count]; ++count; upwd1[k][7] = df1[count]; ++count;
      }
   }

   for (i3=0; i3<n3; ++i3 ) {
      for (i2=0; i2<n2; ++i2 ) {
         k = i2 + n2*i3;
         upwd1[k][1] = df1[count]; ++count; upwd1[k][2] = df1[count]; ++count;
      }
      for (i2=0; i2<n2; ++i2 ) {
         k = i2 + n2*i3;
         upwd1[k][5] = df1[count]; ++count; upwd1[k][6] = df1[count]; ++count;
      }
   }

   //----------------------------------

   count = 0;
   for (i3=0; i3<n3; ++i3 ) {
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i3;
         upwd2[k][0] = df2[count]; ++count; upwd2[k][1] = df2[count]; ++count;
      }
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i3;
         upwd2[k][4] = df2[count]; ++count; upwd2[k][5] = df2[count]; ++count;
      }
   }

   for (i3=0; i3<n3; ++i3 ) {
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i3;
         upwd2[k][3] = df2[count]; ++count; upwd2[k][2] = df2[count]; ++count;
      }
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i3;
         upwd2[k][7] = df2[count]; ++count; upwd2[k][6] = df2[count]; ++count;
      }
   }

   //----------------------------------

   count =0;
   for (i2=0; i2<n2; ++i2 ) {
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i2;
         upwd3[k][0] = df3[count]; ++count; upwd3[k][1] = df3[count]; ++count;
      }
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i2;
         upwd3[k][3] = df3[count]; ++count; upwd3[k][2] = df3[count]; ++count;
      }
   }

   for (i2=0; i2<n2; ++i2 ) {
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i2;
         upwd3[k][4] = df3[count]; ++count; upwd3[k][5] = df3[count]; ++count;
      }
      for (i1=0; i1<n1; ++i1 ) {
         k = i1 + n1*i2;
         upwd3[k][7] = df3[count]; ++count; upwd3[k][6] = df3[count]; ++count;
      }
   }

   

   delete bf1;
   delete bf2;
   delete bf3;
   delete df1;
   delete df2;
   delete df3;

}


//------------------------------------------------------------------------------------------------


void domain_avr_scalar(double c[][8], double cbarx[][8], double cbary[][8], double cbarz[][8],
                       int n1, int n2, int n3, int *mynbr, fstream& out)
{
   // needed by swapbdry()
   double *dptr = NULL;
   int *iptr = NULL;
   int mortar=0, AvgType=1, BdryVarType=2;
   int mynbrSize[12];

   // compute the averages across the edges of the mesh in each subdomain
   avr_x(c, cbarx, n1, n2, n3);
   avr_y(c, cbary, n1, n2, n3);
   avr_z(c, cbarz, n1, n2, n3);

   int nodesX = 2*n1, nodesY = 2*n2, nodesZ = 2*n3; // # of nodes in one direction

   getnbrsize_(mynbr, mynbrSize, nodesX, nodesY, nodesZ);

   double (*bf1) = new double [2*nodesY*nodesZ];
   double (*bf2) = new double [2*nodesX*nodesZ];
   double (*bf3) = new double [2*nodesX*nodesY];
   double (*df1) = new double [2*nodesY*nodesZ];
   double (*df2) = new double [2*nodesX*nodesZ];
   double (*df3) = new double [2*nodesX*nodesY];

   load_buffers(c, bf1, bf2, bf3, n1, n2, n3);
   swapbdry( nodesX, nodesY, nodesZ, mynbr, mynbrSize, bf1, bf2, bf3, df1, df2, df3, BdryVarType, AvgType);
   read_buffers(cbarx, cbary, cbarz, df1, df2, df3, n1, n2, n3);

   delete df3;
   delete df2;
   delete df1;
   delete bf3;
   delete bf2;
   delete bf1;
}


//------------------------------------------------------------------------------------------------


void domain_avr_vector(double qx[][8], double qy[][8], double qz[][8],
                       double qbarx[][8], double qbary[][8], double qbarz[][8],
                       int n1, int n2, int n3, int *mynbr, fstream& out)
{
   double (*qdummy)[8] = new double [n1*n2*n3][8];

   domain_avr_scalar(qx, qbarx, qdummy, qdummy, n1, n2, n3, mynbr, out);
   domain_avr_scalar(qy, qdummy, qbary, qdummy, n1, n2, n3, mynbr, out);
   domain_avr_scalar(qz, qdummy, qdummy, qbarz, n1, n2, n3, mynbr, out);

   delete [] qdummy;
}


//------------------------------------------------------------------------------------------------


void compute_flux( eqn *Ptr, int *mynbr, int &model,
                   double c[][8], double cbarx[][8], double cbary[][8], double cbarz[][8],
                   double qx[][8], double qy[][8], double qz[][8],
                   int& dflag, int& dtype, fstream& out)
{

   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;


   if ( dflag == 0  ) { // NO DIFFUSION
      int elements = n1* n2 * n3;
      for (int k=0; k<elements; ++k) {
	 for (int p=0; p<8; ++p) {
	    qx[k][p] = 0.0; qy[k][p] = 0.0; qz[k][p] = 0.0;
	 }
      }
   }   
   
   else if ( dflag == 1) {

      double (*rhs)[8] = new double [3][8];
      double (*aux)[8] = new double [3][8];
      double tx, ty, tz, coeff, QX, QY, QZ, detx, dety, detz, det;
      int pI, pJ, pindI, pindJ;

      for (int i3=0; i3<n3; ++i3) {
	 for (int i2=0; i2<n2; ++i2) {
	    for (int i1=0; i1<n1; ++i1) {

	       int k = i1 + n1*i2 + n1*n2*i3; // element
	       
	       det = det_JF( Ptr->x1[i1], Ptr->x1[i1+1],
                             Ptr->x2[i2], Ptr->x2[i2+1], 
                             Ptr->x3[i3], Ptr->x3[i3+1] );

	       detx = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
	       dety = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
	       detz = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;


	       for (int p=0; p<8; ++p) {  
		  tx=0.0; ty=0.0; tz=0.0;
		  for (int s=0;s<8; ++s) {
		     tx += c[k][s] * ( Ptr->wwx[s][p] );
		     ty += c[k][s] * ( Ptr->wwy[s][p] );
		     tz += c[k][s] * ( Ptr->wwz[s][p] );
                  }
		  rhs[0][p] = tx * detx; rhs[1][p] = ty * dety; rhs[2][p] = tz * detz;  
	       }

	       
	       // x-sides
	       for (pindJ=0; pindJ<4; ++pindJ) {
		  pJ = Ptr->xside[0][pindJ];
		  tx = 0.0;
		  for (pindI=0; pindI<4; ++pindI) {
		     pI = Ptr->xside[0][pindI];
		     coeff = ( (i1==0) && (mynbr[0]==-1) )  ?  -c[k][pI]  : -cbarx[k][pI];
		     tx += coeff * (Ptr->ww_face[pindI][pindJ]);
		  }
		  rhs[0][pJ] -= tx * detx;
	       }
	       for (pindJ=0; pindJ<4; ++pindJ) {
		  pJ = Ptr->xside[1][pindJ];
		  tx = 0.0;
		  for (pindI=0; pindI<4; ++pindI) {
		     pI = Ptr->xside[1][pindI];
		     coeff = ( (i1==(n1-1)) && (mynbr[3]==-1) )  ?  c[k][pI]  : cbarx[k][pI];
		     tx += coeff * (Ptr->ww_face[pindI][pindJ]);
		  }
		  rhs[0][pJ] -= tx * detx;
	       }

	       // y-sides
	       for (pindJ=0; pindJ<4; ++pindJ) {
		  pJ = Ptr->yside[0][pindJ];
		  ty = 0.0;
		  for (pindI=0; pindI<4; ++pindI) {
		     pI = Ptr->yside[0][pindI];
		     coeff = ( (i2==0) && (mynbr[1]==-1) )  ?  -c[k][pI]  : -cbary[k][pI];
		     ty += coeff * (Ptr->ww_face[pindI][pindJ]);
		  }
		  rhs[1][pJ] -= ty * dety;
	       }
	       for (pindJ=0; pindJ<4; ++pindJ) {
		  pJ = Ptr->yside[1][pindJ];
		  ty = 0.0;
		  for (pindI=0; pindI<4; ++pindI) {
		     pI = Ptr->yside[1][pindI];
		     coeff = ( (i2==(n2-1)) && (mynbr[4]==-1) )  ?  c[k][pI]  : cbary[k][pI];
		     ty += coeff * (Ptr->ww_face[pindI][pindJ]);
		  }
		  rhs[1][pJ] -= ty * dety;
	       }

	       // z-sides
	       for (pindJ=0; pindJ<4; ++pindJ) {
		  pJ = Ptr->zside[0][pindJ];
		  tz = 0.0;
		  for (pindI=0; pindI<4; ++pindI) {
		     pI = Ptr->zside[0][pindI];
		     coeff = ( (i3==0) && (mynbr[2]==-1) )  ?  -c[k][pI]  : -cbarz[k][pI];
		     tz += coeff * (Ptr->ww_face[pindI][pindJ]);
		  }
		  rhs[2][pJ] -= tz * detz;
	       }
	       for (pindJ=0; pindJ<4; ++pindJ) {
		  pJ = Ptr->zside[1][pindJ];
		  tz = 0.0;
		  for (pindI=0; pindI<4; ++pindI) {
		     pI = Ptr->zside[1][pindI];
		     coeff = ( (i3==(n3-1)) && (mynbr[5]==-1) )  ?  c[k][pI]  : cbarz[k][pI];
		     tz += coeff * (Ptr->ww_face[pindI][pindJ]);
		  }
		  rhs[2][pJ] -= tz * detz;
	       }             
	       
	       for (int i=0; i<3; ++i) cholesky_solve(8, Ptr->M, aux[i], rhs[i]);
	       	       
	       for (int p=0; p<8; ++p) {

		  if (dtype==1) {
		     QX = 0.0; QY = 0.0; QZ = 0.0;
		     for (int i=0; i<3; ++i) {
			QX += Ptr->D[0][i] * aux[i][p];
			QY += Ptr->D[1][i] * aux[i][p];
			QZ += Ptr->D[2][i] * aux[i][p];
		     }
		     qx[k][p] = QX/det;
		     qy[k][p] = QY/det;
		     qz[k][p] = QZ/det;
		  }
		  
                  else if (dtype==2) {
		     if (model==1) {
			QX = 0.0; QY = 0.0; QZ = 0.0;
			for (int i=0; i<3; ++i) {
			   QX += Ptr->D_porous[k][0][i] * aux[i][p];
			   QY += Ptr->D_porous[k][1][i] * aux[i][p];
			   QZ += Ptr->D_porous[k][2][i] * aux[i][p];
			}
			qx[k][p] = QX/det;
			qy[k][p] = QY/det;
			qz[k][p] = QZ/det;
		     }
		     if (model==2) {
			qx[k][p] = Ptr->D_fluid * aux[0][p] / det;
			qy[k][p] = Ptr->D_fluid * aux[1][p] / det;
			qz[k][p] = Ptr->D_fluid * aux[2][p] / det;
		     }                     

		  }
                  

	       } // loop over nodes 


	    }
	 }
      }

      delete [] aux;
      delete [] rhs;
   }
   else { cout << "\n dflag = " << dflag << " is not implemented !" << endl << endl; exit(1); }
   

 }


//------------------------------------------------------------------------------------------------


double conc_L2error_sq_subd(int& model, eqn *Ptr, double c[][8], double &t, int& testnumber)
{
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;  
   double norm2 = 0.0, tmp;
   double C[3][3][3], x[3], y[3], z[3];
   double w[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0};

   for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
	 for (int i1=0; i1<n1; ++i1) {
	    int elem  = i1 + n1 * i2  +  n1*n2 * i3;   // element
	    
	    x[0] = Ptr->x1[i1];  x[2] = Ptr->x1[i1+1];  x[1] = ( x[2] + x[0] ) / 2.0;
	    y[0] = Ptr->x2[i2];  y[2] = Ptr->x2[i2+1];  y[1] = ( y[2] + y[0] ) / 2.0;
	    z[0] = Ptr->x3[i3];  z[2] = Ptr->x3[i3+1];  z[1] = ( z[2] + z[0] ) / 2.0;

	    // interpolate to find the values of the concentration at the quadrature points
	    C[0][0][0] = c[elem][0]; 
	    C[2][0][0] = c[elem][1];
	    C[2][2][0] = c[elem][2];
	    C[0][2][0] = c[elem][3];
	    C[0][0][2] = c[elem][4];
	    C[2][0][2] = c[elem][5];
	    C[2][2][2] = c[elem][6];
	    C[0][2][2] = c[elem][7];
	    
	    C[1][0][0] = (c[elem][1] + c[elem][0])/2.0;
	    C[2][1][0] = (c[elem][2] + c[elem][1])/2.0;
	    C[1][2][0] = (c[elem][3] + c[elem][2])/2.0;
	    C[0][1][0] = (c[elem][0] + c[elem][3])/2.0;

	    C[1][0][2] = (c[elem][5] + c[elem][4])/2.0;
	    C[2][1][2] = (c[elem][6] + c[elem][5])/2.0;
	    C[1][2][2] = (c[elem][7] + c[elem][6])/2.0;
	    C[0][1][2] = (c[elem][4] + c[elem][7])/2.0;

	    C[0][0][1] = (c[elem][4] + c[elem][0])/2.0;
	    C[2][0][1] = (c[elem][5] + c[elem][1])/2.0;
	    C[2][2][1] = (c[elem][6] + c[elem][2])/2.0;
	    C[0][2][1] = (c[elem][7] + c[elem][3])/2.0;

	    C[1][0][1] = (C[1][0][2] + C[1][0][0])/2.0;
	    C[2][1][1] = (C[2][1][2] + C[2][1][0])/2.0;
	    C[1][2][1] = (C[1][2][2] + C[1][2][0])/2.0;
	    C[0][1][1] = (C[0][1][2] + C[0][1][0])/2.0;

	    C[1][1][0] = (C[1][2][0] + C[1][0][0])/2.0;
	    C[1][1][2] = (C[1][2][2] + C[1][0][2])/2.0;
	    
	    C[1][1][1] = (C[1][1][2] + C[1][1][0])/2.0;

	    double det = det_JF(Ptr->x1[i1], Ptr->x1[i1+1], Ptr->x2[i2], Ptr->x2[i2+1], Ptr->x3[i3], Ptr->x3[i3+1]);
	    	    
	    for (int i=0; i<3; ++i) {
	       for (int j=0; j<3; ++j) {
		  for (int k=0; k<3; ++k) {
		     tmp = c_true(model, x[i], y[j], z[k], t, testnumber)  -  C[i][j][k];	    
		     norm2 += tmp * tmp * w[i] * w[j] * w[k] * det;
		  }
	       }
	    }


	 }
      }
   }

   return norm2;
}


//------------------------------------------------------------------------------------------------


double conc_L2error_sq(int& model, eqn *Ptr, double c[][8], double &t, int& testnumber)
{
   double norm2 = conc_L2error_sq_subd(model, Ptr, c, t, testnumber);
   gdsumX(6000, &norm2, 1);
   return norm2;
}


double conc_L2error(int& model, eqn *Ptr, double c[][8], double &t, int& testnumber)
{
   double norm2 = conc_L2error_sq_subd(model, Ptr, c, t, testnumber);
   gdsumX(6100, &norm2, 1);
   return sqrt(norm2);
}


//------------------------------------------------------------------------------------------------


double flux_L2error_sq_subd(int& model, eqn *Ptr, double q[][8], int& dir, double &t, int& testnumber)
{
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;  
   double norm2 = 0.0, tmp;
   double Q[3][3][3], x[3], y[3], z[3];
   double w[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0};
   double (*q_true)(int&, eqn*, double, double, double, double, int);

   // pick the function for the true flux corresponding to q
   switch (dir) {
      case 0 : q_true = qx_true; break;
      case 1 : q_true = qy_true; break;
      case 2 : q_true = qz_true; break;
      default: cout << "\nIllegal value for dir in flux_L2error_sq_subd()!  " << dir << endl; exit(1);
   }

   

   for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
	 for (int i1=0; i1<n1; ++i1) {
	    int elem  = i1 + n1 * i2  +  n1*n2 * i3;   // element
	    
	    x[0] = Ptr->x1[i1];  x[2] = Ptr->x1[i1+1];  x[1] = ( x[2] + x[0] ) / 2.0;
	    y[0] = Ptr->x2[i2];  y[2] = Ptr->x2[i2+1];  y[1] = ( y[2] + y[0] ) / 2.0;
	    z[0] = Ptr->x3[i3];  z[2] = Ptr->x3[i3+1];  z[1] = ( z[2] + z[0] ) / 2.0;
	    	    
	    // interpolate to find the values of the flux at the quadrature points
	    Q[0][0][0] = q[elem][0]; 
	    Q[2][0][0] = q[elem][1];
	    Q[2][2][0] = q[elem][2];
	    Q[0][2][0] = q[elem][3];
	    Q[0][0][2] = q[elem][4];
	    Q[2][0][2] = q[elem][5];
	    Q[2][2][2] = q[elem][6];
	    Q[0][2][2] = q[elem][7];
	    
	    Q[1][0][0] = (q[elem][1] + q[elem][0])/2.0;
	    Q[2][1][0] = (q[elem][2] + q[elem][1])/2.0;
	    Q[1][2][0] = (q[elem][3] + q[elem][2])/2.0;
	    Q[0][1][0] = (q[elem][0] + q[elem][3])/2.0;

	    Q[1][0][2] = (q[elem][5] + q[elem][4])/2.0;
	    Q[2][1][2] = (q[elem][6] + q[elem][5])/2.0;
	    Q[1][2][2] = (q[elem][7] + q[elem][6])/2.0;
	    Q[0][1][2] = (q[elem][4] + q[elem][7])/2.0;

	    Q[0][0][1] = (q[elem][4] + q[elem][0])/2.0;
	    Q[2][0][1] = (q[elem][5] + q[elem][1])/2.0;
	    Q[2][2][1] = (q[elem][6] + q[elem][2])/2.0;
	    Q[0][2][1] = (q[elem][7] + q[elem][3])/2.0;

	    Q[1][0][1] = (Q[1][0][2] + Q[1][0][0])/2.0;
	    Q[2][1][1] = (Q[2][1][2] + Q[2][1][0])/2.0;
	    Q[1][2][1] = (Q[1][2][2] + Q[1][2][0])/2.0;
	    Q[0][1][1] = (Q[0][1][2] + Q[0][1][0])/2.0;

	    Q[1][1][0] = (Q[1][2][0] + Q[1][0][0])/2.0;
	    Q[1][1][2] = (Q[1][2][2] + Q[1][0][2])/2.0;
	    
	    Q[1][1][1] = (Q[1][1][2] + Q[1][1][0])/2.0;

	    double det = det_JF(Ptr->x1[i1], Ptr->x1[i1+1], Ptr->x2[i2], Ptr->x2[i2+1], Ptr->x3[i3], Ptr->x3[i3+1]);
	    	    
	    for (int i=0; i<3; ++i) {
	       for (int j=0; j<3; ++j) {
		  for (int k=0; k<3; ++k) {
		     tmp = q_true(model, Ptr, x[i], y[j], z[k], t, testnumber)  -  Q[i][j][k];	    
		     norm2 += tmp * tmp * w[i] * w[j] * w[k] * det;
		  }
	       }
	    }


	 }
      }
   }

   return norm2;
}


//------------------------------------------------------------------------------------------------


double flux_L2error_sq(int& model, eqn *Ptr, double qx[][8], double qy[][8], double qz[][8], double &t, int& testnumber)
{
   int dir;
   double norm2;

   dir = 0; norm2  = flux_L2error_sq_subd(model, Ptr, qx, dir, t, testnumber);
   dir = 1; norm2 += flux_L2error_sq_subd(model, Ptr, qy, dir, t, testnumber);
   dir = 2; norm2 += flux_L2error_sq_subd(model, Ptr, qz, dir, t, testnumber);

   gdsumX(7000, &norm2, 1);
   return norm2;
}


double flux_L2error(int& model, eqn *Ptr, double qx[][8], double qy[][8], double qz[][8], double &t, int& testnumber)
{
   int dir;
   double norm2;

   dir = 0; norm2  = flux_L2error_sq_subd(model, Ptr, qx, dir, t, testnumber);
   dir = 1; norm2 += flux_L2error_sq_subd(model, Ptr, qy, dir, t, testnumber);
   dir = 2; norm2 += flux_L2error_sq_subd(model, Ptr, qz, dir, t, testnumber);

   gdsumX(7100, &norm2, 1);
   return sqrt(norm2);
}


double mass(eqn *Ptr, double c[][8])
{
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;  
   double norm = 0.0, tmp;
   double C[3][3][3];
   double w[3] = {1.0/3.0, 4.0/3.0, 1.0/3.0};

   for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
	 for (int i1=0; i1<n1; ++i1) {
	    int elem  = i1 + n1 * i2  +  n1*n2 * i3;   // element

	    // interpolate to find the values of the concentration at the quadrature points
	    C[0][0][0] = c[elem][0]; 
	    C[2][0][0] = c[elem][1];
	    C[2][2][0] = c[elem][2];
	    C[0][2][0] = c[elem][3];
	    C[0][0][2] = c[elem][4];
	    C[2][0][2] = c[elem][5];
	    C[2][2][2] = c[elem][6];
	    C[0][2][2] = c[elem][7];
	    
	    C[1][0][0] = (c[elem][1] + c[elem][0])/2.0;
	    C[2][1][0] = (c[elem][2] + c[elem][1])/2.0;
	    C[1][2][0] = (c[elem][3] + c[elem][2])/2.0;
	    C[0][1][0] = (c[elem][0] + c[elem][3])/2.0;

	    C[1][0][2] = (c[elem][5] + c[elem][4])/2.0;
	    C[2][1][2] = (c[elem][6] + c[elem][5])/2.0;
	    C[1][2][2] = (c[elem][7] + c[elem][6])/2.0;
	    C[0][1][2] = (c[elem][4] + c[elem][7])/2.0;

	    C[0][0][1] = (c[elem][4] + c[elem][0])/2.0;
	    C[2][0][1] = (c[elem][5] + c[elem][1])/2.0;
	    C[2][2][1] = (c[elem][6] + c[elem][2])/2.0;
	    C[0][2][1] = (c[elem][7] + c[elem][3])/2.0;

	    C[1][0][1] = (C[1][0][2] + C[1][0][0])/2.0;
	    C[2][1][1] = (C[2][1][2] + C[2][1][0])/2.0;
	    C[1][2][1] = (C[1][2][2] + C[1][2][0])/2.0;
	    C[0][1][1] = (C[0][1][2] + C[0][1][0])/2.0;

	    C[1][1][0] = (C[1][2][0] + C[1][0][0])/2.0;
	    C[1][1][2] = (C[1][2][2] + C[1][0][2])/2.0;
	    
	    C[1][1][1] = (C[1][1][2] + C[1][1][0])/2.0;

	    double det = det_JF(Ptr->x1[i1], Ptr->x1[i1+1], Ptr->x2[i2], Ptr->x2[i2+1], Ptr->x3[i3], Ptr->x3[i3+1]);
	    	    
	    for (int i=0; i<3; ++i) {
	       for (int j=0; j<3; ++j) {
		  for (int k=0; k<3; ++k) {
		     tmp = C[i][j][k];	    
		     norm += tmp * w[i] * w[j] * w[k] * det;
		  }
	       }
	    }


	 }
      }
   }

   return norm;


}


//------------------------------------------------------------------------------------------------



void compute_stochparam(int n1, int n2, int n3,
                        double c[][8], double &cweight, fstream& stochIO)
{

   double (*c_exp)[8] = new double [n1*n2*n3][8]; // arrays to store the expectation and the variance
   double (*c_var)[8] = new double [n1*n2*n3][8]; // of the concentration if stochastic loop is used
  

   delete [] c_exp;
   delete [] c_var;

}

// miz17, bag8 : checks to see if stochCount is a member of the RealizationIdx array
bool isRealizationToPlot(int& StochCount, int* RealizationIdx, int* NumOfRealizations) {
   int i;
   for (i=0; i<*NumOfRealizations; i++) {
      // cout << "RealizationIdx[" << i << "]=" << RealizationIdx[i] << endl;
      if (RealizationIdx[i]==StochCount) return true;
   }
   return false;
}

