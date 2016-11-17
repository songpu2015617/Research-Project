#include "LDG.h"

double quad(int& i1, int& i2, int& i3, function wi, function wj, eqn *Ptr)
{
   double integral = quad(wi, wj);   
   integral *= det_JF( Ptr->x1[i1], Ptr->x1[i1+1], Ptr->x2[i2], Ptr->x2[i2+1], Ptr->x3[i3], Ptr->x3[i3+1] );
   return integral;
}

//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------

void LDG_term1(int& i1, int& i2, int& i3, eqn *Ptr,
               double c[][8], double qx[][8], double qy[][8], double qz[][8], double term[8])
{
   int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;
   int n = (Ptr->n1) * (Ptr->n2) * (Ptr->n3);
   double detX = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
   double detY = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
   double detZ = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;
   for (int p=0; p<8; ++p) {
      double termX=0.0, termY=0.0, termZ=0.0;
      for (int m=0; m<8; ++m) {
	 termX += ( c[k][m] * (Ptr->u1[k+m*n]) + qx[k][m] ) * (Ptr->wwx[m][p]);
	 termY += ( c[k][m] * (Ptr->u2[k+m*n]) + qy[k][m] ) * (Ptr->wwy[m][p]);
         termZ += ( c[k][m] * (Ptr->u3[k+m*n]) + qz[k][m] ) * (Ptr->wwz[m][p]);
      }
      termX *= detX; termY *= detY; termZ *= detZ;
      term[p] = termX + termY + termZ;
   }

}

//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------

void LDG_term2(int& model, int& i1, int& i2, int& i3, eqn *Ptr, int *mynbr, 
               double c[][8], double term[8], int& testnumber, double& CurrentTime)
{
    int pI, pJ, pindI, pindJ, x1index, x2index, x3index; 
    double det, u_dot_normal, uvtx_dot_normal, tmp;

    for (pJ=0; pJ<8; ++pJ) term[pJ] = 0.0;

    // ONLY FACES FROM THE GLOBAL DOMAIN BDRY ARE TAKEN

    // check subdomain bdry X=Xmin
    if ( (i1==0) && (mynbr[0]==-1) ) {
       int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;
       int n = (Ptr->n1) * (Ptr->n2) * (Ptr->n3);
       det = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
       u_dot_normal = Ptr->u_normal(k, 0, 3, 7, 4, 'x');
       if ( u_dot_normal<0.0 ) { // inflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->xside[0][pindJ];
	     tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->xside[0][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');
		x1index = i1;
		x2index = ( (pI==0)||(pI==4) ) ? i2 : (i2+1);
		x3index = ( (pI==0)||(pI==3) ) ? i3 : (i3+1);
		tmp += c_inflow(model,
                                Ptr->x1[x1index],
                                Ptr->x2[x2index],
                                Ptr->x3[x3index],
                                CurrentTime,
                                testnumber) * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
	     tmp *= det;
             term[pJ] += tmp;
	  }
       }
       else {                    // outflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->xside[0][pindJ];
	     tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->xside[0][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');
		tmp  += c[k][pI] * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }  
    }
    //------------------------------------------------------
 
    // check subdomain bdry X=Xmax
    if ( (i1==((Ptr->n1)-1)) && (mynbr[3]==-1) ) {
       int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;
       int n = (Ptr->n1) * (Ptr->n2) * (Ptr->n3);
       det = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
       u_dot_normal = Ptr->u_normal(k, 1, 2, 6, 5, 'x');
       if ( u_dot_normal<0.0 ) { // inflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->xside[1][pindJ];
             tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->xside[1][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');
		x1index = i1+1;
		x2index = ( (pI==1)||(pI==5) ) ? i2 : (i2+1);
		x3index = ( (pI==1)||(pI==2) ) ? i3 : (i3+1);
		tmp += c_inflow(model,
                                Ptr->x1[x1index],
                                Ptr->x2[x2index],
                                Ptr->x3[x3index],
                                CurrentTime,
                                testnumber) * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
       else {                    // outflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->xside[1][pindJ];
	     tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->xside[1][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');		
		tmp += c[k][pI] * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
    }
    //------------------------------------------------------
    
    // check subdomain bdry Y=Ymin
    if ( (i2==0) && (mynbr[1]==-1) ) {
       int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;
       int n = (Ptr->n1) * (Ptr->n2) * (Ptr->n3);
       det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
       u_dot_normal = Ptr->u_normal(k, 0, 4, 5, 1, 'y');
       if ( u_dot_normal<0.0 ) { // inflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->yside[0][pindJ];
	     tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->yside[0][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
		x1index = ( (pI==0)||(pI==4) ) ? i1 : (i1+1);
		x2index = i2;
		x3index = ( (pI==0)||(pI==1) ) ? i3 : (i3+1);
		tmp += c_inflow(model,
                                Ptr->x1[x1index],
                                Ptr->x2[x2index],
                                Ptr->x3[x3index],
                                CurrentTime,
                                testnumber) * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
       else {                    // outflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->yside[0][pindJ];
             tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->yside[0][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
		tmp += c[k][pI] * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }  
    }
    //------------------------------------------------------
    

    // check subdomain bdry Y=Ymax
    if ( (i2==((Ptr->n2)-1)) && (mynbr[4]==-1) ) {
       int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;
       int n = (Ptr->n1) * (Ptr->n2) * (Ptr->n3);
       det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
       u_dot_normal = Ptr->u_normal(k, 3, 7, 6, 2, 'y');
       if ( u_dot_normal<0.0 ) { // inflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->yside[1][pindJ];
             tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->yside[1][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
		x1index = ( (pI==3)||(pI==7) ) ? i1 : (i1+1);
		x2index = i2+1;
		x3index = ( (pI==3)||(pI==2) ) ? i3 : (i3+1);
		tmp += c_inflow(model,
                                Ptr->x1[x1index],
                                Ptr->x2[x2index],
                                Ptr->x3[x3index],
                                CurrentTime,
                                testnumber) * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
       else {                    // outflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->yside[1][pindJ];
             tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->yside[1][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
		tmp += c[k][pI] * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
    }
    //------------------------------------------------------
    

    // check subdomain bdry Z=Zmin
    if ( (i3==0) && (mynbr[2]==-1) ) {
       int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;
       int n = (Ptr->n1) * (Ptr->n2) * (Ptr->n3);
       det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;
       u_dot_normal = Ptr->u_normal(k, 0, 1, 2, 3, 'z');
       if ( u_dot_normal<0.0 ) { // inflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->zside[0][pindJ];
	     tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->zside[0][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
		x1index = ( (pI==0)||(pI==3) ) ? i1 : (i1+1);
		x2index = ( (pI==0)||(pI==1) ) ? i2 : (i2+1);
		x3index = i3;
		tmp += c_inflow(model,
                                Ptr->x1[x1index],
                                Ptr->x2[x2index],
                                Ptr->x3[x3index],
                                CurrentTime,
                                testnumber) * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
       else {                    // outflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->zside[0][pindJ];
	     tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->zside[0][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
		tmp += c[k][pI] * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
    }
    //------------------------------------------------------
    

    // check subdomain bdry Z=Zmax
    if ( (i3==((Ptr->n3)-1)) && (mynbr[5]==-1) ) {
       int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;
       int n = (Ptr->n1) * (Ptr->n2) * (Ptr->n3);
       det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;
       u_dot_normal = Ptr->u_normal(k, 4, 5, 6, 7, 'z');
       if ( u_dot_normal<0.0 ) { // inflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->zside[1][pindJ];
	     tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->zside[1][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
		x1index = ( (pI==4)||(pI==7) ) ? i1 : (i1+1);
		x2index = ( (pI==4)||(pI==5) ) ? i2 : (i2+1);
		x3index = i3;
		tmp += c_inflow(model,
                                Ptr->x1[x1index],
                                Ptr->x2[x2index],
                                Ptr->x3[x3index],
                                CurrentTime,
                                testnumber) * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
       else {                    // outflow bdry
	  for (pindJ=0; pindJ<4; ++pindJ) {
	     pJ = Ptr->zside[1][pindJ];
             tmp = 0.0;
	     for (pindI=0; pindI<4; ++pindI) {
		pI = Ptr->zside[1][pindI];
		uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
		tmp += c[k][pI] * uvtx_dot_normal * (Ptr->ww_face[pindI][pindJ]);
	     }
             tmp *= det;
             term[pJ] += tmp;
	  }
       }
    }
    //------------------------------------------------------

}

//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------

void LDG_term3(int& i1, int& i2, int& i3, eqn *Ptr, int *mynbr,
               double qbarx[][8], double qbary[][8], double qbarz[][8],
               double upwd1[][8], double upwd2[][8], double upwd3[][8],
               double c[][8], double term[8])
{
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1 * n2 * n3;
   int k = i1 + (Ptr->n1)*i2 + (Ptr->n1)*(Ptr->n2)*i3;

   int pI, pJ, pindI, pindJ, upwd_el, upwd_vtx;
   double coeff, det, u_dot_normal, uvtx_dot_normal, tmp;

   for (pJ=0; pJ<8; ++pJ) term[pJ] = 0.0;

   // ONLY INTERIOR FACES ARE TAKEN

   // check subdomain bdry X=Xmin
   if (i1>0) {      
      det = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 0, 3, 7, 4, 'x');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->xside[0][pindJ];
         tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->xside[0][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');
	    Ptr->upwd('W', k, upwd_el, pI, upwd_vtx);
	    coeff = (u_dot_normal<0.0) ? c[upwd_el][upwd_vtx] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarx[k][pI] * Ptr->sign[pI][0] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   else if (mynbr[0]>-1) {
      det = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 0, 3, 7, 4, 'x');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->xside[0][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->xside[0][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');
	    coeff = (u_dot_normal<0.0) ? upwd1[i2+n2*i3][pI] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarx[k][pI] * Ptr->sign[pI][0] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   //------------------------------------------------------

   // check subdomain bdry X=Xmax
   if (i1<(n1-1)) {
      det = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 1, 2, 6, 5, 'x');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->xside[1][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->xside[1][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');
	    Ptr->upwd('E', k, upwd_el, pI, upwd_vtx);
	    coeff = (u_dot_normal<0.0) ? c[upwd_el][upwd_vtx] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarx[k][pI] * Ptr->sign[pI][0] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }
   }
   else if (mynbr[3]>-1) {
      det = (Ptr->x2[i2+1] - Ptr->x2[i2]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 1, 2, 6, 5, 'x');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->xside[1][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->xside[1][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'x');
	    coeff = (u_dot_normal<0.0) ? upwd1[i2+n2*i3][pI] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarx[k][pI] * Ptr->sign[pI][0] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }
   }
   //------------------------------------------------------

   // check subdomain bdry Y=Ymin
   if (i2>0) {      
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 0, 4, 5, 1, 'y');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->yside[0][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->yside[0][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
	    Ptr->upwd('S', k, upwd_el, pI, upwd_vtx);
	    coeff = (u_dot_normal<0.0) ? c[upwd_el][upwd_vtx] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbary[k][pI] * Ptr->sign[pI][1] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   else if (mynbr[1]>-1) {
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 0, 4, 5, 1, 'y');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->yside[0][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->yside[0][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
	    coeff = (u_dot_normal<0.0) ? upwd2[i1+n1*i3][pI] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbary[k][pI] * Ptr->sign[pI][1] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   //------------------------------------------------------

    // check subdomain bdry Y=Ymax
   if (i2<(n2-1)) {      
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 3, 7, 6, 2, 'y');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->yside[1][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->yside[1][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
	    Ptr->upwd('N', k, upwd_el, pI, upwd_vtx);
	    coeff = (u_dot_normal<0.0) ? c[upwd_el][upwd_vtx] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbary[k][pI] * Ptr->sign[pI][1] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   else if (mynbr[4]>-1) {
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x3[i3+1] - Ptr->x3[i3]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 3, 7, 6, 2, 'y');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->yside[1][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->yside[1][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'y');
	    coeff = (u_dot_normal<0.0) ? upwd2[i1+n1*i3][pI] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbary[k][pI] * Ptr->sign[pI][1] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   //------------------------------------------------------

   // check subdomain bdry Z=Zmin
   if (i3>0) {      
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 0, 1, 2, 3, 'z');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->zside[0][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->zside[0][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
	    Ptr->upwd('D', k, upwd_el, pI, upwd_vtx);
	    coeff = (u_dot_normal<0.0) ? c[upwd_el][upwd_vtx] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarz[k][pI] * Ptr->sign[pI][2] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   else if (mynbr[2]>-1) {
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 0, 1, 2, 3, 'z');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->zside[0][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->zside[0][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
	    coeff = (u_dot_normal<0.0) ? upwd3[i1+n1*i2][pI] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarz[k][pI] * Ptr->sign[pI][2] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   //------------------------------------------------------

   // check subdomain bdry Z=Zmax
   if (i3<(n3-1)) {      
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 4, 5, 6, 7, 'z');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->zside[1][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->zside[1][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
	    Ptr->upwd('U', k, upwd_el, pI, upwd_vtx);
	    coeff = (u_dot_normal<0.0) ? c[upwd_el][upwd_vtx] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarz[k][pI] * Ptr->sign[pI][2] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   else if (mynbr[5]>-1) {
      det = (Ptr->x1[i1+1] - Ptr->x1[i1]) * (Ptr->x2[i2+1] - Ptr->x2[i2]) / 4.0;
      u_dot_normal = Ptr->u_normal(k, 4, 5, 6, 7, 'z');
      for (pindJ=0; pindJ<4; ++pindJ) {
	 pJ = Ptr->zside[1][pindJ];
	 tmp = 0.0;
	 for (pindI=0; pindI<4; ++pindI) {
	    pI = Ptr->zside[1][pindI];
	    uvtx_dot_normal = Ptr->u_vtx_normal(k, pI, 'z');
	    coeff = (u_dot_normal<0.0) ? upwd3[i1+n1*i2][pI] : c[k][pI];
	    tmp += ( coeff * uvtx_dot_normal + qbarz[k][pI] * Ptr->sign[pI][2] ) * Ptr->ww_face[pindI][pindJ];
	 }
	 tmp *= det;
	 term[pJ] += tmp;
      }      
   }
   //------------------------------------------------------

}


//------------------------------------------------------------------------------------------------

void LDG_term4(int& model, int& i1, int& i2, int& i3, eqn *Ptr, double term[8], int& testnumber, double& CurrentTime)
{
   double det = det_JF(Ptr->x1[i1], Ptr->x1[i1+1], Ptr->x2[i2], Ptr->x2[i2+1], Ptr->x3[i3], Ptr->x3[i3+1]);
   for (int p=0; p<8; ++p) {
      term[p] = 0.0;
      for (int m=0; m<8; ++m) {
	 int x1index = ( (m==0)||(m==3)||(m==7)||(m==4) ) ? i1 : (i1+1);
	 int x2index = ( (m==0)||(m==4)||(m==5)||(m==1) ) ? i2 : (i2+1);
	 int x3index = ( (m==0)||(m==1)||(m==2)||(m==3) ) ? i3 : (i3+1);
	 term[p] +=  f(model, Ptr->x1[x1index], Ptr->x2[x2index], Ptr->x3[x3index], CurrentTime, testnumber, Ptr )
                     *  (Ptr->ww[m][p]);
	 
      }
      term[p] = term[p] * det * (Ptr->porosity);
   }
}

//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
void LDG_rhs(int& model, int& i1, int& i2, int& i3, eqn *Ptr, int *mynbr, double c[][8],
             double upwd1[][8], double upwd2[][8], double upwd3[][8],
             double qx[][8], double qy[][8], double qz[][8],
	     double qbarx[][8], double qbary[][8], double qbarz[][8],
             double TimeStep, double rhs[8], int& testnumber, double& CurrentTime)
{
   int p;
   double *term = new double [8];
   double det = det_JF( Ptr->x1[i1], Ptr->x1[i1+1],
                        Ptr->x2[i2], Ptr->x2[i2+1], 
                        Ptr->x3[i3], Ptr->x3[i3+1] );

   LDG_term1(i1, i2, i3, Ptr, c, qx, qy, qz, rhs);   

   LDG_term2(model, i1, i2, i3, Ptr, mynbr, c, term, testnumber, CurrentTime);
   for (p=0; p<8; ++p) rhs[p] -= term[p];

   LDG_term3(i1, i2, i3, Ptr, mynbr, qbarx, qbary, qbarz, upwd1, upwd2, upwd3, c, term);
   for (p=0; p<8; ++p) rhs[p] -= term[p];

   LDG_term4(model, i1, i2, i3, Ptr, term, testnumber, CurrentTime);
   for (p=0; p<8; ++p) rhs[p] += term[p];

   // rhs is divided by det since 
   // the left hand side of the eqn. (mass matrix) is computed on the reference element
   for (p=0; p<8; ++p) rhs[p] = rhs[p] / det;
  
   delete [] term;
}

//------------------------------------------------------------------------------------------------




//------------------------------------------------------------------------------------------------



void fwd_euler( int& mynod, int& model, int *mynbr, eqn *Ptr, int& limflag, double& alpha, int& dflag, int& dtype,
                double c[][8], double cexp[][8], double cvar[][8], double TimeStep, double TimeMax, int skip,
                fstream& out, int& testnumber, int& ldgerr, int& StochFlag, int& StochCount, double& cweight,
                int *RealizationIdx, int *NumOfRealizations, int *totalloops )
{

   int iterMax = (int) ceil(TimeMax / TimeStep);
   int i1, i2, i3, p,j, k, framecheck=-1, current_frame_number=0;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1*n2*n3; // number of elements

   double c_err_L2, c_err_inf, q_err_L2=0.0, q_err_inf=0.0, tmp, sqrt_tmp, CurrentTime=0.0, quad_weight;
   double InitialMass, FinalMass;

   double *M = new double [64];
   double *RHS = new double[8]; 
   double *c_update = new double[8];

   // array to store the old time step
   double (*cold)[8] = new double [n][8];

   double (*qx)[8] = new double [n1*n2*n3][8]; // x-component of the flux
   double (*qy)[8] = new double [n1*n2*n3][8]; // y-component of the flux
   double (*qz)[8] = new double [n1*n2*n3][8]; // z-component of the flux

   double (*barx)[8] = new double [n1*n2*n3][8]; // averages along x-
   double (*bary)[8] = new double [n1*n2*n3][8]; // averages along y-
   double (*barz)[8] = new double [n1*n2*n3][8]; // averages along z-

   // arrays to store the concentration values of the neighbouring processors along the interface
   double (*upwd1)[8] = new double [n2*n3][8];
   double (*upwd2)[8] = new double [n1*n3][8];
   double (*upwd3)[8] = new double [n1*n2][8]; 


   fstream solnOUT, stochIO;
   /*
      solnOUT is used to store the solution if StochFlag=0
              or a particular realization if StochFlag=1

      stochIO is used to store the mean and the variance (StochFlag=1)
   */      

   // create the initial frame  (t = t0)
   char initialframe[40];
   char initialstochframe[40];
   if (StochFlag==0) {
      strcpy(initialframe,"transport_output/t_out_00_0000.plt");
      frame_name(mynod, 0, initialframe);
      solnOUT.open( initialframe, ios::out );
      print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
   }
   else {
      if ( isRealizationToPlot(StochCount,RealizationIdx,NumOfRealizations) ) {
         strcpy(initialframe,"transport_output/t_mc0000_00_0000.plt");
         realization_frame_name(mynod, StochCount, 0, initialframe);
         solnOUT.open( initialframe, ios::out );
         print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
      }
      // miz17, bag8 - fill cexp and cvar then write to file
      bool first_realization;
      strcpy(initialstochframe,"transport_output/stoch_00_0000.plt");
      frame_name(mynod, 0, initialstochframe);
      if (StochCount==0) {
         first_realization = true;
         read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
      }
      else {
         first_realization = false;
         stochIO.open( initialstochframe, ios::in | ios::binary );
         read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
         // if final realization, cvar = cvar - cexp^2
         if (*totalloops==(StochCount+1)) 
            subtract_square_exp(Ptr, cexp, cvar);
         stochIO.close();
      }

      if (*totalloops==(StochCount+1)) {
         stochIO.open( initialstochframe, ios::out );
         print_stochastic(Ptr, cexp, cvar, stochIO, current_frame_number, CurrentTime, mynod);
      }
      else {
         stochIO.open( initialstochframe, ios::out | ios::binary );
         write_binary( Ptr, cexp, cvar, stochIO);
      }
   }
   if ( solnOUT.is_open() ) solnOUT.close();
   if ( stochIO.is_open() ) stochIO.close();
   ++current_frame_number;
   ++framecheck;

   //for(int i=0; i<8; ++i) for(int j=0; j<8; ++j) M[i+8*j] = Ptr->ww[i][j];
   /*** use Cholesky factorization for the mass matrix ***/ 
   //char uplo='U';
   //int dim=8, nrhs=1, info;
   //dposv_(uplo, dim, nrhs, M, dim, RHS, dim, info);
   /******************************************************/

   /*
       Forward Euler scheme

       y'(t) = f(t,y)

       K = f( t[n-1] , y[n-1] )
       y[n] = y[n-1] + dt * K

   */

   domain_avr_scalar(c, barx, bary, barz, n1, n2, n3, mynbr, out);       
   compute_flux(Ptr, mynbr, model, c, barx, bary, barz, qx, qy, qz, dflag, dtype, out);          
   domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
   interface_val(c, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

   // compute the contributions to the errors at CurrentTime = 0.0
   if (ldgerr==1) {
      quad_weight = 0.5;
      tmp = conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber);	 
      c_err_L2 = tmp * quad_weight;
      c_err_inf = sqrt(tmp);
      if (dflag!=0) {
	 tmp = flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber);
	 q_err_L2 = tmp * quad_weight;
	 q_err_inf = sqrt(tmp);
      }
   }

   InitialMass = mass(Ptr,c);

   for (int iter=1; iter<=iterMax; ++iter) { // time loop

      // bag8 : debug
      if (mynod==0) cout << "fwd_euler iter = " << iter << endl;

      for (k=0; k<n; ++k) for (p=0; p<8; ++p) cold[k][p] = c[k][p];
      
      // dhv      
      // print_output(Ptr, barx, bary, barz, out);

      for (i3=0; i3<n3; ++i3) {
         for (i2=0; i2<n2; ++i2) {
            for (i1=0; i1<n1; ++i1) {

               k  = i1 + n1 * i2  +  n1*n2 * i3;   // element

	       LDG_rhs(model, i1, i2, i3, Ptr, mynbr, cold,
	               upwd1, upwd2, upwd3, qx, qy, qz, barx, bary, barz,
                       TimeStep, RHS, testnumber, CurrentTime);

               cholesky_solve(8, Ptr->M, c[k], RHS);   // find K

            }
         }
      }

      for (k=0; k<n; ++k) for (p=0; p<8; ++p) c[k][p] = cold[k][p] + TimeStep * c[k][p];
      if (limflag==1) limiter(alpha, n1, n2, n3, mynbr, c, out);      
      
      domain_avr_scalar(c, barx, bary, barz, n1, n2, n3, mynbr, out);       
      compute_flux(Ptr, mynbr, model, c, barx, bary, barz, qx, qy, qz, dflag, dtype, out);          
      domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
      interface_val(c, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

      CurrentTime += TimeStep;

      // compute the contributions to the errors at CurrentTime = iter * TimeStep
      if (ldgerr==1) {
	 quad_weight = (iter==iterMax)  ?  0.5  :  1.0;
	 tmp = conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber);	 
	 c_err_L2 += tmp * quad_weight;
	 sqrt_tmp = sqrt(tmp); c_err_inf = ( c_err_inf > sqrt_tmp ) ? c_err_inf : sqrt_tmp;
         if (dflag!=0) {
	    tmp = flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber);
	    q_err_L2 += tmp * quad_weight;
	    sqrt_tmp = sqrt(tmp); q_err_inf = ( q_err_inf > sqrt_tmp ) ? q_err_inf : sqrt_tmp;
         }
      }

      // create a frame for t > t0
      if (framecheck == skip)
      {
         char frame[40];
         char stochframe[40];
         if (StochFlag==0) {
            strcpy(frame,"transport_output/t_out_00_0000.plt");
            frame_name(mynod, current_frame_number, frame);
            solnOUT.open( frame, ios::out );
            print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
         }
         else {
            if ( isRealizationToPlot(StochCount,RealizationIdx,NumOfRealizations) ) {
               strcpy(frame,"transport_output/t_mc0000_00_0000.plt");
               realization_frame_name(mynod, StochCount, current_frame_number, frame);
               solnOUT.open( frame, ios::out );
               print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
            }
            // miz17, bag8 - fill cexp and cvar then write to file
            bool first_realization;
            strcpy(stochframe,"transport_output/stoch_00_0000.plt");
            frame_name(mynod, current_frame_number, stochframe);
            if (StochCount==0) {
               first_realization = true;
               read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
            }
            else {
               first_realization = false;
               stochIO.open(stochframe, ios::in | ios::binary);
               read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
               // if final realization, cvar = cvar - cexp^2
               if (*totalloops==(StochCount+1)) 
                  subtract_square_exp(Ptr, cexp, cvar);
               stochIO.close();
            }

            if (*totalloops==(StochCount+1)) {
               stochIO.open( stochframe, ios::out );
               print_stochastic(Ptr, cexp, cvar, stochIO, current_frame_number, CurrentTime, mynod);
            }
            else {
               stochIO.open( stochframe, ios::out | ios::binary );
               write_binary( Ptr, cexp, cvar, stochIO);
            }
         }
         if ( solnOUT.is_open() ) solnOUT.close();
         if ( stochIO.is_open() ) stochIO.close();
         ++current_frame_number;
         framecheck = -1;
      }
      ++framecheck;

   } // iter

   if (ldgerr==1) {
      c_err_L2  =  sqrt( c_err_L2 * TimeMax / iterMax );
      q_err_L2  =  (dflag!=0) ? sqrt( q_err_L2 * TimeMax / iterMax ) : q_err_L2;
   }
   if (ldgerr==2) {
      tmp =  conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber); c_err_L2 = sqrt(tmp);
      if (dflag!=0) {
	 tmp =  flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber); q_err_L2 = sqrt(tmp);
      }
   }

   FinalMass = mass(Ptr,c);

   //  print the errors
   if (mynod==0) {
      if (ldgerr==1) {
	 cout << "LDG errors:" << endl; 
	 cout << "C L2(0,T;L2)   C Linf(0,T;L2)   Q L2(0,T;L2)   Q Linf(0,T;L2)" << endl;
	 cout << scientific << c_err_L2 << "   " << c_err_inf << "     " 
              << q_err_L2 << "   " << q_err_inf << endl << endl;
	 //cout << "L2(0,T;L2) error of the concentration   = " << scientific << c_err_L2 << endl;
	 //cout << "L2(0,T;L2) error of the flux            = " << scientific << q_err_L2 << endl;
	 //cout << "Linf(0,T;L2) error of the concentration = " << scientific << c_err_inf << endl;
	 //cout << "Linf(0,T;L2) error of the flux          = " << scientific << q_err_inf << endl;
	 //cout << endl;
      }
      if (ldgerr==2) {
	 cout << "LDG errors in L2-norm at time t = " << CurrentTime << endl;
	 cout << " C               Q" << endl;
	 cout << scientific << c_err_L2 << "   " <<  q_err_L2 << endl << endl;
      }
      cout << "Mass balance: " << scientific << InitialMass << " " << FinalMass << endl << endl;
   }
   

   delete [] upwd3;
   delete [] upwd2;
   delete [] upwd1;
   delete [] barz;
   delete [] bary;
   delete [] barx;
   delete [] qz;
   delete [] qy;
   delete [] qx;
   delete [] cold;
   delete [] c_update;
   delete [] RHS;
   delete [] M;
}




//------------------------------------------------------------------------------------------------




void rk2( int& mynod, int& model, int *mynbr, eqn *Ptr, int& limflag, double& alpha, int& dflag, int& dtype,
          double c[][8], double cexp[][8], double cvar[][8], double TimeStep, double TimeMax, int skip,
	      fstream& out, int& testnumber, int& ldgerr,int& StochFlag, int& StochCount, double& cweight,
          int *RealizationIdx, int *NumOfRealizations, int *totalloops )
{

   int iterMax = (int) ceil(TimeMax / TimeStep);
   int i1, i2, i3, p,j, k, framecheck=-1, current_frame_number=0;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1*n2*n3; // number of elements

   double c_err_L2, c_err_inf, q_err_L2=0.0, q_err_inf=0.0, tmp, sqrt_tmp, CurrentTime=0.0, quad_weight;
   double InitialMass, FinalMass;

   double *M = new double [64];
   double *RHS = new double [8]; 
   double *c_update = new double [8];

   // array needed for the Runge-Kutta method
   double (*rk_stage)[8] = new double [n][8];

   // array to store the old time step
   double (*cold)[8] = new double [n][8];

   double (*qx)[8] = new double [n1*n2*n3][8]; // x-component of the flux
   double (*qy)[8] = new double [n1*n2*n3][8]; // y-component of the flux
   double (*qz)[8] = new double [n1*n2*n3][8]; // z-component of the flux

   double (*barx)[8] = new double [n1*n2*n3][8]; // averages along x-
   double (*bary)[8] = new double [n1*n2*n3][8]; // averages along y-
   double (*barz)[8] = new double [n1*n2*n3][8]; // averages along z-

   // arrays to store the concentration values of the neighbouring processors along the interface
   double (*upwd1)[8] = new double [n2*n3][8];
   double (*upwd2)[8] = new double [n1*n3][8];
   double (*upwd3)[8] = new double [n1*n2][8]; 

   fstream solnOUT, stochIO;
   /*
      solnOUT is used to store the solution if StochFlag=0
              or a particular realization if StochFlag=1

      stochIO is used to store the mean and the variance (StochFlag=1)
   */ 
   // create the initial frame  (t = t0)
   char initialframe[40];
   char initialstochframe[40];
   if (StochFlag==0) {
      strcpy(initialframe,"transport_output/t_out_00_0000.plt");
      frame_name(mynod, 0, initialframe);
      solnOUT.open( initialframe, ios::out );
      print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
   }
   else {
      if ( isRealizationToPlot(StochCount,RealizationIdx,NumOfRealizations) ) {
         strcpy(initialframe,"transport_output/t_mc0000_00_0000.plt");
         realization_frame_name(mynod, StochCount, 0, initialframe);
         solnOUT.open( initialframe, ios::out );
         print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
      }
      // miz17, bag8 - fill cexp and cvar then write to file
      bool first_realization;
      strcpy(initialstochframe,"transport_output/stoch_00_0000.plt");
      frame_name(mynod, 0, initialstochframe);
      if (StochCount==0) {
         first_realization = true;
         read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
      }
      else {
         first_realization = false;
         stochIO.open( initialstochframe, ios::in | ios::binary);
         read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
         // if final realization, cvar = cvar - cexp^2
         if (*totalloops==(StochCount+1)) 
            subtract_square_exp(Ptr, cexp, cvar);
         stochIO.close();
      }
      if (*totalloops==(StochCount+1)) {
         stochIO.open( initialstochframe, ios::out );
         print_stochastic(Ptr, cexp, cvar, stochIO, current_frame_number, CurrentTime, mynod);
      }
      else {
         stochIO.open( initialstochframe, ios::out | ios::binary );
         write_binary( Ptr, cexp, cvar, stochIO);
      }
   }
   if ( solnOUT.is_open() ) solnOUT.close();
   if ( stochIO.is_open() ) stochIO.close();
   ++current_frame_number;
   ++framecheck;


   //for(int i=0; i<8; ++i) for(int j=0; j<8; ++j) M[i+8*j] = Ptr->ww[i][j];
   /*** use Cholesky factorization for the mass matrix ***/ 
   //char uplo='U';
   //int dim=8, nrhs=1, info;
   //dposv_(uplo, dim, nrhs, M, dim, RHS, dim, info);
   /******************************************************/

   /*
       RK2 scheme

       y'(t) = f(t,y)

       K1 = f( t[n-1] , y[n-1] )
       K2 = f( t[n-1] + dt/2 , y[n-1] + (dt/2) * K1 )

       y[n] = y[n-1] + dt * K2

   */
   
   domain_avr_scalar(c, barx, bary, barz, n1, n2, n3, mynbr, out);
   compute_flux(Ptr, mynbr, model, c, barx, bary, barz, qx, qy, qz, dflag, dtype, out);
   domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
   interface_val(c, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

   
   // compute the contributions to the errors at CurrentTime = 0.0
   if (ldgerr==1) {
      quad_weight = 0.5;
      tmp = conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber);	 
      c_err_L2 = tmp * quad_weight;
      c_err_inf = sqrt(tmp);
      if (dflag!=0) {
	 tmp = flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber);
	 q_err_L2 = tmp * quad_weight;
	 q_err_inf = sqrt(tmp);
      }
   }

   InitialMass = mass(Ptr, c);


   for (int iter=1; iter<=iterMax; ++iter) { // time loop

      // bag8 : debug
      if (mynod==0) cout << "rk2 iter = " << iter << endl;

      for (k=0; k<n; ++k) for (p=0; p<8; ++p) cold[k][p] = c[k][p];

      //------------------------------ 1st stage ---------------------------------      
	 
      for (i3=0; i3<n3; ++i3) {
	 for (i2=0; i2<n2; ++i2) {
	    for (i1=0; i1<n1; ++i1) {
	       
	       k  = i1 + n1 * i2  +  n1*n2 * i3;   // element
	       
	       LDG_rhs(model, i1, i2, i3, Ptr, mynbr, cold,
		       upwd1, upwd2, upwd3, qx, qy, qz, barx, bary, barz,
		       TimeStep, RHS, testnumber, CurrentTime);
	       
	       cholesky_solve(8, Ptr->M, c[k], RHS);        // find K1 

	    }
	 }
      }
      //--------------------------------------------------------------------------

      double rk_t = CurrentTime + 0.5 * TimeStep;
      for (k=0; k<n; ++k) for (p=0; p<8; ++p) rk_stage[k][p] = cold[k][p] + 0.5 * TimeStep * c[k][p];
            
      domain_avr_scalar(rk_stage, barx, bary, barz, n1, n2, n3, mynbr, out);
      compute_flux(Ptr, mynbr, model, rk_stage, barx, bary, barz, qx, qy, qz, dflag, dtype, out);      
      domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
      interface_val(rk_stage, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

      //------------------------------ 2nd stage ---------------------------------
      for (i3=0; i3<n3; ++i3) {
	 for (i2=0; i2<n2; ++i2) {
	    for (i1=0; i1<n1; ++i1) {

	       k  = i1 + n1 * i2  +  n1*n2 * i3;   // element

	       LDG_rhs(model, i1, i2, i3, Ptr, mynbr, rk_stage,
		       upwd1, upwd2, upwd3, qx, qy, qz, barx, bary, barz,
		       TimeStep, RHS, testnumber, rk_t);
	       
	       cholesky_solve(8, Ptr->M, c[k], RHS);        // find K2

	    }
	 }
      }
      //--------------------------------------------------------------------------

      for (k=0; k<n; ++k) for (p=0; p<8; ++p) c[k][p] = cold[k][p] + TimeStep * c[k][p];
      if (limflag==1) limiter(alpha, n1, n2, n3, mynbr, c, out);

      domain_avr_scalar(c, barx, bary, barz, n1, n2, n3, mynbr, out);       
      compute_flux(Ptr, mynbr, model, c, barx, bary, barz, qx, qy, qz, dflag, dtype, out);
      domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
      interface_val(c, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

      CurrentTime += TimeStep;

      // compute the contributions to the errors at CurrentTime = iter * TimeStep
      if (ldgerr==1) {
	 quad_weight = (iter==iterMax)  ?  0.5  :  1.0;
	 tmp = conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber);	 
	 c_err_L2 += tmp * quad_weight;
	 sqrt_tmp = sqrt(tmp); c_err_inf = ( c_err_inf > sqrt_tmp ) ? c_err_inf : sqrt_tmp;
	 if (dflag!=0) {
	    tmp = flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber);
	    q_err_L2 += tmp * quad_weight;
	    sqrt_tmp = sqrt(tmp); q_err_inf = ( q_err_inf > sqrt_tmp ) ? q_err_inf : sqrt_tmp;
         }
      }


      // create a frame for t > t0
      if (framecheck == skip)
      {
         char frame[40];
         char stochframe[40];
         if (StochFlag==0) {
            strcpy(frame,"transport_output/t_out_00_0000.plt");
            frame_name(mynod, current_frame_number, frame);
            solnOUT.open( frame, ios::out );
            print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
         }
         else {
            if ( isRealizationToPlot(StochCount,RealizationIdx,NumOfRealizations) ) {
               strcpy(frame,"transport_output/t_mc0000_00_0000.plt");
               realization_frame_name(mynod, StochCount, current_frame_number, frame);
               solnOUT.open( frame, ios::out );
               print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
            }
            // miz17, bag8 - fill cexp and cvar then write to file
            bool first_realization;
            strcpy(stochframe,"transport_output/stoch_00_0000.plt");
            frame_name(mynod, current_frame_number, stochframe);
            if (StochCount==0) {
               first_realization = true;
               read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
            }
            else {
               first_realization = false;
               stochIO.open(stochframe, ios::in | ios::binary);
               read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
               // if final realization, cvar = cvar - cexp^2
               if (*totalloops==(StochCount+1)) 
                  subtract_square_exp(Ptr, cexp, cvar);
               stochIO.close();
            }
            if (*totalloops==(StochCount+1)) {
               stochIO.open( stochframe, ios::out );
               print_stochastic(Ptr, cexp, cvar, stochIO, current_frame_number, CurrentTime, mynod);
            }
            else {
               stochIO.open( stochframe, ios::out | ios::binary );
               write_binary( Ptr, cexp, cvar, stochIO);
            }
         }
         if ( solnOUT.is_open() ) solnOUT.close();
         if ( stochIO.is_open() ) stochIO.close();
         ++current_frame_number;
         framecheck = -1;
      }
      ++framecheck;

   } // iter

   FinalMass = mass(Ptr, c);

   if (ldgerr==1) {
      c_err_L2  =  sqrt( c_err_L2 * TimeMax / iterMax );
      q_err_L2  =  (dflag!=0) ? sqrt( q_err_L2 * TimeMax / iterMax ) : q_err_L2;
   }
   if (ldgerr==2) {
      tmp =  conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber); c_err_L2 = sqrt(tmp);
      if (dflag!=0) {
	 tmp =  flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber); q_err_L2 = sqrt(tmp);
      }
   }

   //  print the errors
   if (mynod==0) {
      if (ldgerr==1) {
	 cout << "LDG errors:" << endl; 
	 cout << "C L2(0,T;L2)   C Linf(0,T;L2)   Q L2(0,T;L2)   Q Linf(0,T;L2)" << endl;
	 cout << scientific << c_err_L2 << "   " << c_err_inf << "     " 
              << q_err_L2 << "   " << q_err_inf << endl << endl;
	 //cout << "L2(0,T;L2) error of the concentration   = " << scientific << c_err_L2 << endl;
	 //cout << "L2(0,T;L2) error of the flux            = " << scientific << q_err_L2 << endl;
	 //cout << "Linf(0,T;L2) error of the concentration = " << scientific << c_err_inf << endl;
	 //cout << "Linf(0,T;L2) error of the flux          = " << scientific << q_err_inf << endl;
	 //cout << endl;
      }
      if (ldgerr==2) {
	 cout << "LDG errors in L2-norm at time t = " << CurrentTime << endl;
	 cout << " C               Q" << endl;
	 cout << scientific << c_err_L2 << "   " <<  q_err_L2 << endl << endl;
      }
      cout << "Mass balance: " << scientific << InitialMass << " " << FinalMass << endl << endl;
   }

   delete [] upwd3;
   delete [] upwd2;
   delete [] upwd1;
   delete [] barz;
   delete [] bary;
   delete [] barx;
   delete [] qz;
   delete [] qy;
   delete [] qx;
   delete [] cold;
   delete [] rk_stage;
   delete [] c_update;
   delete [] RHS;
   delete [] M;

}




//------------------------------------------------------------------------------------------------




void rk4( int& mynod, int& model, int *mynbr, eqn *Ptr, int& limflag, double& alpha, int& dflag, int& dtype,
          double c[][8], double cexp[][8], double cvar[][8], double TimeStep, double TimeMax, int skip,
	      fstream& out, int& testnumber, int& ldgerr, int& StochFlag, int& StochCount, double& cweight,
          int *RealizationIdx, int *NumOfRealizations, int *totalloops )
{

   int iterMax = (int) ceil(TimeMax / TimeStep);
   int i1, i2, i3, p,j, k, framecheck=-1, current_frame_number=0;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1*n2*n3; // number of elements

   double c_err_L2, c_err_inf, q_err_L2=0.0, q_err_inf=0.0, tmp, sqrt_tmp, CurrentTime=0.0, quad_weight;
   double InitialMass, FinalMass;

   double *M = new double [64];
   double *RHS = new double [8]; 
   double *c_update = new double [8];

   // arrays needed for the Runge-Kutta method
   double (*rk_stage)[8] = new double [n][8];
   double (*rk_slope)[8] = new double [n][8];

   // array to store the old time step
   double (*cold)[8] = new double [n][8];

   double (*qx)[8] = new double [n1*n2*n3][8]; // x-component of the flux
   double (*qy)[8] = new double [n1*n2*n3][8]; // y-component of the flux
   double (*qz)[8] = new double [n1*n2*n3][8]; // z-component of the flux

   double (*barx)[8] = new double [n1*n2*n3][8]; // averages along x-
   double (*bary)[8] = new double [n1*n2*n3][8]; // averages along y-
   double (*barz)[8] = new double [n1*n2*n3][8]; // averages along z-

   // arrays to store the concentration values of the neighbouring processors along the interface
   double (*upwd1)[8] = new double [n2*n3][8];
   double (*upwd2)[8] = new double [n1*n3][8];
   double (*upwd3)[8] = new double [n1*n2][8]; 

   fstream solnOUT, stochIO;
   /*
      solnOUT is used to store the solution if StochFlag=0
              or a particular realization if StochFlag=1

      stochIO is used to store the mean and the variance (StochFlag=1)
   */ 
   // create the initial frame  (t = t0)
   char initialframe[40];
   char initialstochframe[40];
   if (StochFlag==0) {
      strcpy(initialframe,"transport_output/t_out_00_0000.plt");
      frame_name(mynod, 0, initialframe);
      solnOUT.open( initialframe, ios::out );
      print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
   }
   else {
      if ( isRealizationToPlot(StochCount,RealizationIdx,NumOfRealizations) ) {
         strcpy(initialframe,"transport_output/t_mc0000_00_0000.plt");
         realization_frame_name(mynod, StochCount, 0, initialframe);
         solnOUT.open( initialframe, ios::out );
         print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
      }
      // miz17, bag8 - fill cexp and cvar then write to file
      bool first_realization;
      strcpy(initialstochframe,"transport_output/stoch_00_0000.plt");
      frame_name(mynod, 0, initialstochframe);
      if (StochCount==0) {
         first_realization = true;
         read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
      }
      else {
         first_realization = false;
         stochIO.open(initialstochframe, ios::in | ios::binary);
         read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
         // if final realization, cvar = cvar - cexp^2
         if (*totalloops==(StochCount+1)) 
            subtract_square_exp(Ptr, cexp, cvar);
         stochIO.close();
      }
      if (*totalloops==(StochCount+1)) {
         stochIO.open( initialstochframe, ios::out );
         print_stochastic(Ptr, cexp, cvar, stochIO, current_frame_number, CurrentTime, mynod);
      }
      else {
         stochIO.open( initialstochframe, ios::out | ios::binary );
         write_binary( Ptr, cexp, cvar, stochIO);
      }
   }
   if ( solnOUT.is_open() ) solnOUT.close();
   if ( stochIO.is_open() ) stochIO.close();
   ++current_frame_number;
   ++framecheck;


   //for(int i=0; i<8; ++i) for(int j=0; j<8; ++j) M[i+8*j] = Ptr->ww[i][j];
   /*** use Cholesky factorization for the mass matrix ***/ 
   //char uplo='U';
   //int dim=8, nrhs=1, info;
   //dposv_(uplo, dim, nrhs, M, dim, RHS, dim, info);
   /******************************************************/

   /*
       RK4 scheme

       y'(t) = f(t,y)

       K1 = f( t[n-1] , y[n-1] )
       K2 = f( t[n-1] + dt/2 , y[n-1] + (dt/2) * K1 )
       K3 = f( t[n-1] + dt/2 , y[n-1] + (dt/2) * K2 )
       K4 = f( t[n-1] +  dt  , y[n-1] +  dt * K3 )

       y[n] = y[n-1] + dt * ( K1/6 + K2/3 + K3/3 + K4/6 ) 

   */

   domain_avr_scalar(c, barx, bary, barz, n1, n2, n3, mynbr, out);       
   compute_flux(Ptr, mynbr, model, c, barx, bary, barz, qx, qy, qz, dflag, dtype, out);
   domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
   interface_val(c, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

   // compute the contributions to the errors at CurrentTime = 0.0
   if (ldgerr==1) {
      quad_weight = 0.5;
      tmp = conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber);	 
      c_err_L2 = tmp * quad_weight;
      c_err_inf = sqrt(tmp);
      if (dflag!=0) {
	 tmp = flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber);
	 q_err_L2 = tmp * quad_weight;
	 q_err_inf = sqrt(tmp);
      }
   }

   InitialMass = mass(Ptr,c);

   for (int iter=1; iter<=iterMax; ++iter) { // time loop

      // bag8 : debug
      if (mynod==0) cout << "rk4 iter = " << iter << endl;

      for (k=0; k<n; ++k) for (p=0; p<8; ++p) cold[k][p] = c[k][p];

      //------------------------------ 1st stage ---------------------------------      	 
      for (i3=0; i3<n3; ++i3) {
	 for (i2=0; i2<n2; ++i2) {
	    for (i1=0; i1<n1; ++i1) {
	       
	       k  = i1 + n1 * i2  +  n1*n2 * i3;   // element
	       
	       LDG_rhs(model, i1, i2, i3, Ptr, mynbr, cold,
		       upwd1, upwd2, upwd3, qx, qy, qz, barx, bary, barz,
		       TimeStep, RHS, testnumber, CurrentTime);
	       
	       cholesky_solve(8, Ptr->M, rk_slope[k], RHS);        // find K1 

	    }
	 }
      }
      //--------------------------------------------------------------------------

      double rk_t = CurrentTime + 0.5 * TimeStep;
      for (k=0; k<n; ++k) for (p=0; p<8; ++p) {
	 c[k][p] = rk_slope[k][p] / 6.0;
	 rk_stage[k][p] = cold[k][p] + 0.5 * TimeStep * rk_slope[k][p];
      }
      domain_avr_scalar(rk_stage, barx, bary, barz, n1, n2, n3, mynbr, out);       
      compute_flux(Ptr, mynbr, model, rk_stage, barx, bary, barz, qx, qy, qz, dflag, dtype, out);
      domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
      interface_val(rk_stage, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

      //------------------------------ 2nd stage ---------------------------------
      for (i3=0; i3<n3; ++i3) {
	 for (i2=0; i2<n2; ++i2) {
	    for (i1=0; i1<n1; ++i1) {

	       k  = i1 + n1 * i2  +  n1*n2 * i3;   // element

	       LDG_rhs(model, i1, i2, i3, Ptr, mynbr, rk_stage,
		       upwd1, upwd2, upwd3, qx, qy, qz, barx, bary, barz,
		       TimeStep, RHS, testnumber, rk_t);
	       
	       cholesky_solve(8, Ptr->M, rk_slope[k], RHS);        // find K2

	    }
	 }
      }
      //--------------------------------------------------------------------------

      for (k=0; k<n; ++k) for (p=0; p<8; ++p) {
	 c[k][p] += rk_slope[k][p] / 3.0;
	 rk_stage[k][p] = cold[k][p] + 0.5 * TimeStep * rk_slope[k][p];
      }
      domain_avr_scalar(rk_stage, barx, bary, barz, n1, n2, n3, mynbr, out);       
      compute_flux(Ptr, mynbr, model, rk_stage, barx, bary, barz, qx, qy, qz, dflag, dtype, out);      
      domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
      interface_val(rk_stage, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

      //------------------------------ 3rd stage ---------------------------------
      for (i3=0; i3<n3; ++i3) {
	 for (i2=0; i2<n2; ++i2) {
	    for (i1=0; i1<n1; ++i1) {

	       k  = i1 + n1 * i2  +  n1*n2 * i3;   // element

	       LDG_rhs(model, i1, i2, i3, Ptr, mynbr, rk_stage,
		       upwd1, upwd2, upwd3, qx, qy, qz, barx, bary, barz,
		       TimeStep, RHS, testnumber, rk_t);
	       
	       cholesky_solve(8, Ptr->M, rk_slope[k], RHS);        // find K3

	    }
	 }
      }
      //--------------------------------------------------------------------------

      rk_t = CurrentTime + TimeStep;
      for (k=0; k<n; ++k) for (p=0; p<8; ++p) {
	 c[k][p] += rk_slope[k][p] / 3.0;
	 rk_stage[k][p] = cold[k][p] + TimeStep * rk_slope[k][p];
      }
      domain_avr_scalar(rk_stage, barx, bary, barz, n1, n2, n3, mynbr, out);       
      compute_flux(Ptr, mynbr, model, rk_stage, barx, bary, barz, qx, qy, qz, dflag, dtype, out);      
      domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
      interface_val(rk_stage, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

      //------------------------------ 4th stage ---------------------------------
      for (i3=0; i3<n3; ++i3) {
	 for (i2=0; i2<n2; ++i2) {
	    for (i1=0; i1<n1; ++i1) {

	       k  = i1 + n1 * i2  +  n1*n2 * i3;   // element

	       LDG_rhs(model, i1, i2, i3, Ptr, mynbr, rk_stage,
		       upwd1, upwd2, upwd3, qx, qy, qz, barx, bary, barz,
		       TimeStep, RHS, testnumber, rk_t);
	       
	       cholesky_solve(8, Ptr->M, rk_slope[k], RHS);        // find K4

	    }
	 }
      }
      //--------------------------------------------------------------------------

      for (k=0; k<n; ++k) for (p=0; p<8; ++p) {
	 c[k][p] += rk_slope[k][p] / 6.0;
         c[k][p] *= TimeStep;
	 c[k][p] += cold[k][p];
      }
      if (limflag==1) limiter(alpha, n1, n2, n3, mynbr, c, out);

      domain_avr_scalar(c, barx, bary, barz, n1, n2, n3, mynbr, out);       
      compute_flux(Ptr, mynbr, model, c, barx, bary, barz, qx, qy, qz, dflag, dtype, out);
      domain_avr_vector(qx, qy, qz, barx, bary, barz, n1, n2, n3, mynbr, out);
      interface_val(c, upwd1, upwd2, upwd3, n1, n2, n3, mynbr, out);

      CurrentTime += TimeStep;

      // compute the contributions to the errors at CurrentTime = iter * TimeStep
      if (ldgerr==1) {
	 quad_weight = (iter==iterMax)  ?  0.5  :  1.0;
	 tmp = conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber);	 
	 c_err_L2 += tmp * quad_weight;
	 sqrt_tmp = sqrt(tmp); c_err_inf = ( c_err_inf > sqrt_tmp ) ? c_err_inf : sqrt_tmp;
	 if (dflag!=0) {
	    tmp = flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber);
	    q_err_L2 += tmp * quad_weight;
	    sqrt_tmp = sqrt(tmp); q_err_inf = ( q_err_inf > sqrt_tmp ) ? q_err_inf : sqrt_tmp;
         }
      }
      // create a frame for t > t0
      if (framecheck == skip)
      {
         char frame[40];
         char stochframe[40];
         if (StochFlag==0) {
            strcpy(frame,"transport_output/t_out_00_0000.plt");
            frame_name(mynod, current_frame_number, frame);
            solnOUT.open( frame, ios::out );
            print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
         }
         else {
            if ( isRealizationToPlot(StochCount,RealizationIdx,NumOfRealizations) ) {
               strcpy(frame,"transport_output/t_mc0000_00_0000.plt");
               realization_frame_name(mynod, StochCount, current_frame_number, frame);
               solnOUT.open( frame, ios::out );
               print_output(model, Ptr, c, solnOUT, current_frame_number, CurrentTime, testnumber, mynod);
            }
            // miz17, bag8 - fill cexp and cvar then write to file
            bool first_realization;
            strcpy(stochframe,"transport_output/stoch_00_0000.plt");
            frame_name(mynod, current_frame_number, stochframe);
            if (StochCount==0) {
               first_realization = true;
               read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
            }
            else {
               first_realization = false;
               stochIO.open(stochframe, ios::in | ios::binary);
               read_add_stochastic(Ptr, c, cexp, cvar, cweight, stochIO, first_realization, mynod);
               // if final realization, cvar = cvar - cexp^2
               if (*totalloops==(StochCount+1)) 
                  subtract_square_exp(Ptr, cexp, cvar);
               stochIO.close();
            }
            if (*totalloops==(StochCount+1)) {
               stochIO.open( stochframe, ios::out );
               print_stochastic(Ptr, cexp, cvar, stochIO, current_frame_number, CurrentTime, mynod);
            }
            else {
               stochIO.open( stochframe, ios::out | ios::binary );
               write_binary( Ptr, cexp, cvar, stochIO);
            }
         }
         if ( solnOUT.is_open() ) solnOUT.close();
         if ( stochIO.is_open() ) stochIO.close();
         ++current_frame_number;
         framecheck = -1;
      }
      ++framecheck;      

   } // iter

   if (ldgerr==1) {
      c_err_L2  =  sqrt( c_err_L2 * TimeMax / iterMax );
      q_err_L2  =  (dflag!=0) ? sqrt( q_err_L2 * TimeMax / iterMax ) : q_err_L2;
   }
   if (ldgerr==2) {
      tmp =  conc_L2error_sq(model, Ptr, c, CurrentTime, testnumber); c_err_L2 = sqrt(tmp);
      if (dflag!=0) {
	 tmp =  flux_L2error_sq(model, Ptr, qx, qy, qz, CurrentTime, testnumber); q_err_L2 = sqrt(tmp);
      }
   }

   FinalMass = mass(Ptr,c);

   //  print the errors
   if (mynod==0) {
      if (ldgerr==1) {
	 cout << "LDG errors:" << endl; 
	 cout << "C L2(0,T;L2)   C Linf(0,T;L2)   Q L2(0,T;L2)   Q Linf(0,T;L2)" << endl;
	 cout << scientific << c_err_L2 << "   " << c_err_inf << "     " 
              << q_err_L2 << "   " << q_err_inf << endl << endl;
	 //cout << "L2(0,T;L2) error of the concentration   = " << scientific << c_err_L2 << endl;
	 //cout << "L2(0,T;L2) error of the flux            = " << scientific << q_err_L2 << endl;
	 //cout << "Linf(0,T;L2) error of the concentration = " << scientific << c_err_inf << endl;
	 //cout << "Linf(0,T;L2) error of the flux          = " << scientific << q_err_inf << endl;
	 //cout << endl;
      }
      if (ldgerr==2) {
	 cout << "LDG errors in L2-norm at time t = " << CurrentTime << endl;
	 cout << " C               Q" << endl;
	 cout << scientific << c_err_L2 << "   " <<  q_err_L2 << endl << endl;
      }
      cout << "Mass balance: " << scientific << InitialMass << " " << FinalMass << endl << endl;
   }



   delete [] upwd3;
   delete [] upwd2;
   delete [] upwd1;
   delete [] barz;
   delete [] bary;
   delete [] barx;
   delete [] qz;
   delete [] qy;
   delete [] qx;
   delete [] cold;
   delete [] rk_slope;
   delete [] rk_stage;
   delete [] c_update;
   delete [] RHS;
   delete [] M;

}




//------------------------------------------------------------------------------------------------

