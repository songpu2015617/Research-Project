#ifndef _EQN_H
#define _EQN_H

#include "common.h"
#include "basis.h"
#include "quad.h"

class eqn {

 public:

   int n1, n2, n3; // number of partitions along the coord. axes
   double *x1, *x2, *x3; // mesh points
   double *u1, *u2, *u3; // velocity values
   double bdry[3][2];    // computational domain boundary

   double D[3][3];  // diffusion tensor

   // for more 'realistic' diffusion
   double (*D_porous)[3][3]; 
   double D_fluid;

   // makes sense for the porous medium only
   double porosity;

   //indices of the nodes on a given side
   int xside[2][4], yside[2][4], zside[2][4];

   // used for psi psi' terms
   int nbr_vtx[8][3];
   double sign[8][3];

   double M[64];
   double ww[8][8];
   double wwx[8][8];
   double wwy[8][8];
   double wwz[8][8];
   double ww_face[4][4];
   
   eqn( int nx, int ny, int nz, double domainbdry[3][2],
	double *crd1, double *crd2, double *crd3,
	double *vel1, double *vel2, double *vel3 );

   void init();
   double volume(int i1, int i2, int i3);
   double areaX(int i2, int i3);
   double areaY(int i1, int i3);
   double areaZ(int i1, int i2);
   double integral_psipsi_cell(int i1, int i2, int i3);
   double integral_psipsi_side(int i, int j, char side);
   double integral_psipsi_side(int i, int j, int side);
   double integral_psidpsi_cell(int i1, int i2, int i3, int vtx, int crd);
   void upwd(char dir, int el, int& nbr_el, int& ind1, int& ind2, int& ind3, int& ind4);
   void upwd(char dir, int el, int& nbr_el, int vtx, int& nbr_vtx);
   int bdry_type(int element, int vtx1, int vtx2, int vtx3, int vtx4, char side);
   double u_normal(int element, int vtx1, int vtx2, int vtx3, int vtx4, char side);
   double u_vtx_normal(int element, int vtx, char side);

};


#endif
