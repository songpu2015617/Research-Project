#include "eqn.h"

eqn::eqn( int nx, int ny, int nz, double domainbdry[3][2],
              double *crd1, double *crd2, double *crd3,
              double *vel1, double *vel2, double *vel3 )
{
   n1=nx; n2=ny; n3=nz;
   x1=crd1; x2=crd2; x3=crd3;
   u1=vel1; u2=vel2; u3=vel3;

   for (int i=0; i<3; ++i) for (int j=0; j<2; ++j) bdry[i][j] = domainbdry[i][j];

   for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) D[i][j] = (i==j) ? 1.0 : 0.0;

   porosity = 1.0;

   init();
}

void eqn::init()
{
   xside[0][0] = 0; xside[0][1] = 3; xside[0][2] = 7; xside[0][3] = 4;
   xside[1][0] = 1; xside[1][1] = 2; xside[1][2] = 6; xside[1][3] = 5;
   yside[0][0] = 0; yside[0][1] = 4; yside[0][2] = 5; yside[0][3] = 1;
   yside[1][0] = 3; yside[1][1] = 7; yside[1][2] = 6; yside[1][3] = 2;
   zside[0][0] = 0; zside[0][1] = 1; zside[0][2] = 2; zside[0][3] = 3;
   zside[1][0] = 4; zside[1][1] = 5; zside[1][2] = 6; zside[1][3] = 7;


   //     x nbr              y nbr              z nbr
   nbr_vtx[0][0] = 1; nbr_vtx[0][1] = 3; nbr_vtx[0][2] = 4;
   nbr_vtx[1][0] = 0; nbr_vtx[1][1] = 2; nbr_vtx[1][2] = 5;
   nbr_vtx[2][0] = 3; nbr_vtx[2][1] = 1; nbr_vtx[2][2] = 6;
   nbr_vtx[3][0] = 2; nbr_vtx[3][1] = 0; nbr_vtx[3][2] = 7;
   nbr_vtx[4][0] = 5; nbr_vtx[4][1] = 7; nbr_vtx[4][2] = 0;
   nbr_vtx[5][0] = 4; nbr_vtx[5][1] = 6; nbr_vtx[5][2] = 1;
   nbr_vtx[6][0] = 7; nbr_vtx[6][1] = 5; nbr_vtx[6][2] = 2;
   nbr_vtx[7][0] = 6; nbr_vtx[7][1] = 4; nbr_vtx[7][2] = 3;

   //  x                  y                  z
   sign[0][0] = -1.0; sign[0][1] = -1.0; sign[0][2] = -1.0;
   sign[1][0] =  1.0; sign[1][1] = -1.0; sign[1][2] = -1.0;
   sign[2][0] =  1.0; sign[2][1] =  1.0; sign[2][2] = -1.0;
   sign[3][0] = -1.0; sign[3][1] =  1.0; sign[3][2] = -1.0;
   sign[4][0] = -1.0; sign[4][1] = -1.0; sign[4][2] =  1.0;
   sign[5][0] =  1.0; sign[5][1] = -1.0; sign[5][2] =  1.0;
   sign[6][0] =  1.0; sign[6][1] =  1.0; sign[6][2] =  1.0;
   sign[7][0] = -1.0; sign[7][1] =  1.0; sign[7][2] =  1.0;


   function w[8], dw[8];

   w[0] = w0; w[1] = w1; w[2] = w2; w[3] = w3;
   w[4] = w4; w[5] = w5; w[6] = w6; w[7] = w7;
   for (int i=0; i<8; ++i) for (int j=0; j<8; ++j) { ww[i][j] = quad(w[i], w[j]); M[i+8*j] = ww[i][j]; }

   dw[0] = w0x; dw[1] = w1x; dw[2] = w2x; dw[3] = w3x;
   dw[4] = w4x; dw[5] = w5x; dw[6] = w6x; dw[7] = w7x;
   for (int i=0; i<8; ++i) for (int j=0; j<8; ++j) wwx[i][j] = quad(w[i], dw[j]);

   dw[0] = w0y; dw[1] = w1y; dw[2] = w2y; dw[3] = w3y;
   dw[4] = w4y; dw[5] = w5y; dw[6] = w6y; dw[7] = w7y;
   for (int i=0; i<8; ++i) for (int j=0; j<8; ++j) wwy[i][j] = quad(w[i], dw[j]);
   
   dw[0] = w0z; dw[1] = w1z; dw[2] = w2z; dw[3] = w3z;
   dw[4] = w4z; dw[5] = w5z; dw[6] = w6z; dw[7] = w7z;
   for (int i=0; i<8; ++i) for (int j=0; j<8; ++j) wwz[i][j] = quad(w[i], dw[j]);

   for (int i=0; i<4; ++i) for (int j=0; j<4; ++j) ww_face[i][j] = quad(2, -1.0, w[i], w[j]);
   //ww_face[0][0] = 4.0/9.0; ww_face[0][1] = 2.0/9.0; ww_face[0][2] = 1.0/9.0; ww_face[0][3] = 2.0/9.0;
   //ww_face[1][0] = 2.0/9.0; ww_face[1][1] = 4.0/9.0; ww_face[1][2] = 2.0/9.0; ww_face[1][3] = 1.0/9.0;
   //ww_face[2][0] = 1.0/9.0; ww_face[2][1] = 2.0/9.0; ww_face[2][2] = 4.0/9.0; ww_face[2][3] = 2.0/9.0;
   //ww_face[3][0] = 2.0/9.0; ww_face[3][1] = 1.0/9.0; ww_face[3][2] = 2.0/9.0; ww_face[3][3] = 4.0/9.0;
}

double eqn::volume(int i1, int i2, int i3)
{ return (x1[i1+1] - x1[i1]) * (x2[i2+1] - x2[i2]) * (x3[i3+1] - x3[i3]); }

double eqn::areaX(int i2, int i3) { return (x2[i2+1] - x2[i2]) * (x3[i3+1] - x3[i3]); }

double eqn::areaY(int i1, int i3) { return (x1[i1+1] - x1[i1]) * (x3[i3+1] - x3[i3]); }

double eqn::areaZ(int i1, int i2) { return (x1[i1+1] - x1[i1]) * (x2[i2+1] - x2[i2]); }

double eqn::integral_psipsi_cell(int i1, int i2, int i3) { return volume(i1, i2, i3)/8.0; }

double eqn::integral_psipsi_side(int i, int j, char side)
{
   double tmp;
   if (side=='x') { tmp = areaX(i,j)/4.0; }
   if (side=='y') { tmp = areaY(i,j)/4.0; }
   if (side=='z') { tmp = areaZ(i,j)/4.0; }
   return tmp;
}

double eqn::integral_psipsi_side(int i, int j, int side)
{
   double tmp;
   if (side==0) { tmp = areaX(i,j)/4.0; }
   if (side==1) { tmp = areaY(i,j)/4.0; }
   if (side==2) { tmp = areaZ(i,j)/4.0; }
   return tmp;
}

double eqn::integral_psidpsi_cell(int i1, int i2, int i3, int vtx, int crd)
{  double val;
   if (crd==0) val = sign[vtx][crd] * areaX(i2, i3) / 4.0;
   if (crd==1) val = sign[vtx][crd] * areaY(i1, i3) / 4.0;
   if (crd==2) val = sign[vtx][crd] * areaZ(i1, i2) / 4.0;
   return val;
}

void eqn::upwd(char dir, int el, int& nbr_el, int& ind1, int& ind2, int& ind3, int& ind4)
{

   if (dir=='S') { nbr_el=el-n1; ind1=3; ind2=2; ind3=7; ind4=6; }
   if (dir=='N') { nbr_el=el+n1; ind1=0; ind2=1; ind3=4; ind4=5; }

   if (dir=='W') { nbr_el=el-1; ind1=1; ind2=2; ind3=5; ind4=6; }
   if (dir=='E') { nbr_el=el+1; ind1=0; ind2=3; ind3=4; ind4=7; }

   if (dir=='D') { nbr_el=el-n1*n2; ind1=4; ind2=5; ind3=7; ind4=6; }
   if (dir=='U') { nbr_el=el+n1*n2; ind1=0; ind2=1; ind3=3; ind4=2; }
}

void eqn::upwd(char dir, int el, int& nbr_el, int vtx, int& nbr_vtx)
{
   if (dir=='S') { nbr_el=el-n1;
      switch (vtx) {
         case 0: nbr_vtx=3; break;
         case 1: nbr_vtx=2; break;
         case 4: nbr_vtx=7; break;
         case 5: nbr_vtx=6; break;
      }
   }
   if (dir=='N') { nbr_el=el+n1;
      switch (vtx) {
         case 3: nbr_vtx=0; break;
         case 2: nbr_vtx=1; break;
         case 7: nbr_vtx=4; break;
         case 6: nbr_vtx=5; break;
      }
   }
   if (dir=='W') { nbr_el=el-1;
      switch (vtx) {
         case 0: nbr_vtx=1; break;
         case 3: nbr_vtx=2; break;
         case 4: nbr_vtx=5; break;
         case 7: nbr_vtx=6; break;
      }
   }
   if (dir=='E') { nbr_el=el+1;
      switch (vtx) {
         case 1: nbr_vtx=0; break;
         case 2: nbr_vtx=3; break;
         case 5: nbr_vtx=4; break;
         case 6: nbr_vtx=7; break;
      }
   }
   if (dir=='D') { nbr_el=el-n1*n2;
      switch (vtx) {
         case 0: nbr_vtx=4; break;
         case 1: nbr_vtx=5; break;
         case 3: nbr_vtx=7; break;
         case 2: nbr_vtx=6; break;
      }
   }
   if (dir=='U') { nbr_el=el+n1*n2;
      switch (vtx) {
         case 4: nbr_vtx=0; break;
         case 5: nbr_vtx=1; break;
         case 7: nbr_vtx=3; break;
         case 6: nbr_vtx=2; break;
      }
   }
}

int eqn::bdry_type(int element, int vtx1, int vtx2, int vtx3, int vtx4, char side)
{
   int n = n1*n2*n3; // total number of elements
   double avr;
   int bdry;

   if (side=='x') {
      double normalx = (vtx1==0) ? -1.0 : 1.0;
      avr = normalx * ( u1[element + vtx1 * n] + u1[element + vtx2 * n] +
                        u1[element + vtx3 * n] + u1[element + vtx4 * n] );
   }
   if (side=='y') {
      double normaly = (vtx1==0) ? -1.0 : 1.0;
      avr = normaly * ( u2[element + vtx1 * n] + u2[element + vtx2 * n] +
                        u2[element + vtx3 * n] + u2[element + vtx4 * n] );
   }
   if (side=='z') {
      double normalz = (vtx1==0) ? -1.0 : 1.0;
      avr = normalz * ( u3[element + vtx1 * n] + u3[element + vtx2 * n] +
                        u3[element + vtx3 * n] + u3[element + vtx4 * n] );
   }

   bdry = (avr<0.0) ? -1 : 1;
   return bdry;
}

double eqn::u_normal(int element, int vtx1, int vtx2, int vtx3, int vtx4, char side)
{
   int n = n1*n2*n3; // total number of elements
   double avr;

   if (side=='x') {
      double normalx = (vtx1==0) ? -1.0 : 1.0;
      avr = normalx * ( u1[element + vtx1 * n] + u1[element + vtx2 * n] +
                        u1[element + vtx3 * n] + u1[element + vtx4 * n] ) / 4.0;
   }
   if (side=='y') {
      double normaly = (vtx1==0) ? -1.0 : 1.0;
      avr = normaly * ( u2[element + vtx1 * n] + u2[element + vtx2 * n] +
                        u2[element + vtx3 * n] + u2[element + vtx4 * n] ) / 4.0;
   }
   if (side=='z') {
      double normalz = (vtx1==0) ? -1.0 : 1.0;
      avr = normalz * ( u3[element + vtx1 * n] + u3[element + vtx2 * n] +
                        u3[element + vtx3 * n] + u3[element + vtx4 * n] ) / 4.0;
   }
   return avr;
}

double eqn::u_vtx_normal(int element, int vtx, char side)
{
   int n = n1*n2*n3; // total number of elements
   double proj;
   if (side=='x') {
      double normalx = ((vtx==0)||(vtx==3)||(vtx==7)||(vtx==4)) ? -1.0 : 1.0;
      proj = normalx * u1[element + vtx*n];
   }
   if (side=='y') {
      double normaly = ((vtx==0)||(vtx==4)||(vtx==5)||(vtx==1)) ? -1.0 : 1.0;
      proj = normaly * u2[element + vtx*n];
   }
   if (side=='z') {
      double normalz = ((vtx==0)||(vtx==1)||(vtx==2)||(vtx==3)) ? -1.0 : 1.0;
      proj = normalz * u3[element + vtx*n];
   }
   return proj;
}

