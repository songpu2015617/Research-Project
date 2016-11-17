#include "limiter.h"

double minmod(double a1, double a2, double a3)
{
   int s1 = sign(a1), s2 = sign(a2), s3 = sign(a3);
   double res;
   if ( (s1 == s2) && (s2 == s3) ) {
      double val1 = fabs(a1), val2 = fabs(a2), val3 = fabs(a3);
      res = (val1 < val2) ? val1 : val2;
      res = (res < val3) ? res : val3;
      res *= s1;
   }
   else { res = 0.0; }   
   return res;   
}



//------------------------------------------------------------------------------------------------


double cell_avr(double c[8])
{
   double avr = 0.0;
   for (int i=0; i<8; ++i) avr += c[i];
   avr = avr / 8.0;
   return avr;
}


//------------------------------------------------------------------------------------------------


void cell_avr_array(int& n1, int& n2, int& n3, double c[][8], double avr_c[])
{
   for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
	 for (int i1=0; i1<n1; ++i1) {
	    int k = i1 + n1*i2 + n1*n2*i3;
            avr_c[k] = cell_avr( c[k] );
	 }
      }
   }
}


//------------------------------------------------------------------------------------------------


void limiter(double& alpha, int& n1, int& n2, int& n3, int *mynbr, double c[][8], fstream& out)
  
{
   // recover the values at the nodes C from the face averages F
   //                    C = M * F / 3
   double M[8][6] =
      { 2.0, -1.0,  2.0, -1.0,  2.0, -1.0,
       -1.0,  2.0,  2.0, -1.0,  2.0, -1.0,
       -1.0,  2.0, -1.0,  2.0,  2.0, -1.0,
	2.0, -1.0, -1.0,  2.0,  2.0, -1.0,
	2.0, -1.0,  2.0, -1.0, -1.0,  2.0,
       -1.0,  2.0,  2.0, -1.0, -1.0,  2.0,
       -1.0,  2.0, -1.0,  2.0, -1.0,  2.0,
	2.0, -1.0, -1.0,  2.0, -1.0,  2.0 };

   double *avr_c =  new double [n1*n2*n3];
   cell_avr_array(n1, n2, n3, c, avr_c);

   // print the cell averages
   /*
   for (int i=0; i<(n1*n2*n3); ++i) {
      cout << c[i][0] << " " << c[i][1] << " " << c[i][2] << " " << c[i][3] << " "
	   << c[i][4] << " " << c[i][5] << " " << c[i][6] << " " << c[i][7] << "   "
	   << avr_c[i] << endl;
   }
   cout << "-------------" << endl;
   */
   
   double (*bf1) = new double [2*n2*n3];
   double (*bf2) = new double [2*n1*n3];
   double (*bf3) = new double [2*n1*n2];
   double (*df1) = new double [2*n2*n3];
   double (*df2) = new double [2*n1*n3];
   double (*df3) = new double [2*n1*n2];

   int count, i1, i2, i3;   

   count = 0;
   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         i1=0;    bf1[count]       = cell_avr(c[i1 + n1*i2 + n1*n2*i3]);
         i1=n1-1; bf1[count+n2*n3] = cell_avr(c[i1 + n1*i2 + n1*n2*i3]);
         ++count;
      }
   }

   count = 0;
   for (i3=0; i3<n3; ++i3) {
      for (i1=0; i1<n1; ++i1) {
	 i2=0;    bf2[count]       = cell_avr(c[i1 + n1*i2 + n1*n2*i3]);
         i2=n2-1; bf2[count+n1*n3] = cell_avr(c[i1 + n1*i2 + n1*n2*i3]);
         ++count;
      }
   }

   count = 0;
   for (i2=0; i2<n2; ++i2) {
      for (i1=0; i1<n1; ++i1) {
	 i3=0;    bf3[count]       = cell_avr(c[i1 + n1*i2 + n1*n2*i3]);
         i3=n3-1; bf3[count+n1*n2] = cell_avr(c[i1 + n1*i2 + n1*n2*i3]);
         ++count;
      }
   }

   // needed by swapbdry()
   double *dptr = NULL;
   int *iptr = NULL;
   int mortar=0, AvgType=4, BdryVarType=2;
   int mynbrSize[12];
   // -----------------------


   //if (mynbr[3]==3) {for (int i=0; i<2*n2*n3; ++i) bf1[i]=1.0; }
   //for (int i=0; i<2*n2*n3; ++i) out << bf1[i] << endl;
   //out << "******" << endl;
   //for (int i=0; i<2*n2*n3; ++i) out << df1[i] << endl;
   //out << "******" << endl;
   //out << "******" << endl;
   //out << "******" << endl;

   getnbrsize_(mynbr, mynbrSize, n1, n2, n3);
   swapbdry( n1, n2, n3, mynbr, mynbrSize, bf1, bf2, bf3, df1, df2, df3, BdryVarType, AvgType);

   //for (int i=0; i<2*n2*n3; ++i) out << bf1[i] << endl;
   //out << "******" << endl;
   //for (int i=0; i<2*n2*n3; ++i) out << df1[i] << endl;
   //out << "******" << endl;

   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
	 for (i1=0; i1<n1; ++i1) {

	    int k = i1 + n1*i2 + n1*n2*i3;

            double a1, a2, f[6];            
            f[0] = arithmetic_mean(c[k][0], c[k][3], c[k][4], c[k][7]);
            f[1] = arithmetic_mean(c[k][1], c[k][2], c[k][5], c[k][6]);
            f[2] = arithmetic_mean(c[k][0], c[k][1], c[k][4], c[k][5]);
            f[3] = arithmetic_mean(c[k][3], c[k][2], c[k][7], c[k][6]);
            f[4] = arithmetic_mean(c[k][0], c[k][1], c[k][3], c[k][2]);
            f[5] = arithmetic_mean(c[k][4], c[k][5], c[k][7], c[k][6]);

	    
            //--------------------------------------------------------------------------------------
            a1 = ( i1 > 0 )      ?    alpha*(avr_c[k] - avr_c[k-1]) : 
                                      (  ( mynbr[0]>-1 ) ? alpha*(avr_c[k] - df1[i2+n2*i3]) : 0.0  );
            a2 = ( i1 < (n1-1) ) ?    alpha*(avr_c[k+1] - avr_c[k]) : 
                                      (  ( mynbr[3]>-1 ) ? alpha*(df1[i2+n2*i3+n2*n3] - avr_c[k]) : 0.0  );
            f[0] = avr_c[k] - minmod( avr_c[k]-f[0], a1, a2);
            f[1] = avr_c[k] + minmod( f[1]-avr_c[k], a1, a2);
            //--------------------------------------------------------------------------------------           
            a1 = ( i2 > 0 )      ?    alpha*(avr_c[k] - avr_c[k-n1]) :
	                              (  ( mynbr[1]>-1 ) ? alpha*(avr_c[k] - df2[i1+n1*i3]) : 0.0  );
            a2 = ( i2 < (n2-1) ) ?    alpha*(avr_c[k+n1]-avr_c[k])   :
	                              (  ( mynbr[4]>-1 ) ? alpha*(df2[i1+n1*i3+n1*n3] - avr_c[k]) : 0.0  );
            f[2] = avr_c[k] - minmod( avr_c[k]-f[2], a1, a2);
            f[3] = avr_c[k] + minmod( f[3]-avr_c[k], a1, a2);
            //--------------------------------------------------------------------------------------
            a1 = ( i3 > 0 )      ?    alpha*(avr_c[k] - avr_c[k-n1*n2*i3]) :
	                              (  ( mynbr[2]>-1 ) ? alpha*(avr_c[k] - df3[i1+n1*i2]) : 0.0  );
            a2 = ( i3 < (n3-1) ) ?    alpha*(avr_c[k+n1*n2]-avr_c[k])   :
	                              (  ( mynbr[5]>-1 ) ? alpha*(df3[i1+n1*i2+n1*n2] - avr_c[k]) : 0.0  );
            f[4] = avr_c[k] - minmod(avr_c[k]-f[4], a1, a2);
            f[5] = avr_c[k] + minmod(f[5]-avr_c[k], a1, a2);
            //--------------------------------------------------------------------------------------
	    

            
	    //double sum_before = 0.0;
	    //for (int i=0; i<8; ++i) sum_before += c[k][i];
            

            
            for (int i=0; i<8; ++i) {
	       double tmp = 0.0;
               for (int j=0; j<6; ++j)   tmp += M[i][j] * f[j];
               c[k][i] = tmp/3.0;
            }

            /*
	    // recover c[k][i] from the faces sharing the vertex i
            // this approach is very diffusive
	    c[k][0] = ( f[0] + f[2] + f[4] ) / 3.0;
	    c[k][1] = ( f[1] + f[2] + f[4] ) / 3.0;
	    c[k][2] = ( f[1] + f[3] + f[4] ) / 3.0;
	    c[k][3] = ( f[0] + f[3] + f[4] ) / 3.0;
	    c[k][4] = ( f[0] + f[2] + f[5] ) / 3.0;
	    c[k][5] = ( f[1] + f[2] + f[5] ) / 3.0;
	    c[k][6] = ( f[1] + f[3] + f[5] ) / 3.0;
	    c[k][7] = ( f[0] + f[3] + f[5] ) / 3.0;
            */
            
	    //double sum_after = 0.0;
	    //for (int i=0; i<8; ++i) sum_after += c[k][i];
	    //cout << k << " " << sum_before << " " << sum_after << endl;
            
	 }
      }
   }
   //cout << "--------" << endl;

   delete df3;
   delete df2;
   delete df1;
   delete bf3;
   delete bf2;
   delete bf1;
   delete avr_c;
}

