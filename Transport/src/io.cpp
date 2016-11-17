#include "io.h"


//------------------------------------------------------------------------------------------------


void read_file(ifstream& in, double bdry[3][2],
               int& tflag, int& velflag, int& testnumber,
               int& limflag, double& alpha, int& ldgerr, int& dflag, int& dtype, double D[],
               double& d_molecular, double& d_longitudinal, double& d_transverse,
               double& porosity, double& d_visc,
               int& IntegrationMethod, double& T, double& dt, int& skip)

{
   string ln;

   //------------------- Go to the "Domain Size" section in the input file -------------------
   while ( in.good() ) {
      getline(in,ln);
      size_t found = ln.find("x1_min");
      if ( found != string::npos ) break;
   }
   size_t number_begin = 0, number_end = ln.find_first_of(" ");
   string number = ln.substr(number_begin, number_end);
   bdry[0][0] = strtod( number.c_str(), NULL ); // store x1_min
   in >> bdry[0][1]; getline(in,ln);            // store x1_max and go to the next line
 
   for (int i=1; i<3; ++i) { 
      for (int j=0; j<2; ++j) {
	 in >> bdry[i][j]; getline(in,ln);
      }
   }
   //-----------------------------------------------------------------------------------------

   //---------------------------------- Go to Transport --------------------------------------
   while ( in.good() ) { getline(in,ln); if (ln[0] == ':') break; }   

   in >> tflag; getline(in,ln);
   in >> velflag; getline(in,ln);
   in >> testnumber; getline(in,ln);
   in >> IntegrationMethod; getline(in,ln);
   in >> T; getline(in,ln);
   in >> dt; getline(in,ln);
   in >> skip; getline(in,ln);
   in >> dflag; getline(in,ln);
   in >> dtype; getline(in,ln);
   for (int i=0; i<9; ++i) { in >> D[i]; getline(in,ln); }
   in >> d_molecular; getline(in,ln);
   in >> d_longitudinal; getline(in,ln);
   in >> d_transverse; getline(in,ln);
   in >> porosity; getline(in,ln);
   in >> d_visc; getline(in,ln);
   in >> limflag; getline(in,ln);
   in >> alpha; getline(in,ln);
   in >> ldgerr; getline(in,ln);
}
//-----------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------


void readtparam_(int& tflag, int& velflag, int& testnumber, int& dtype)
{

   string ln;

   ifstream in("input");

   // Go to Transport
   while ( in.good() ) { getline(in,ln); if (ln[0] == ':') break; }
   in >> tflag; getline(in,ln);
   in >> velflag; getline(in,ln);
   in >> testnumber; getline(in,ln);
   getline(in,ln); getline(in,ln); getline(in,ln); getline(in,ln); getline(in,ln);
   in >> dtype; getline(in,ln);

   in.close();

}


//------------------------------------------------------------------------------------------------


void read_file_dbg(ifstream& in, double vel[])
{
   string ln;

   // Go to DEBUG DATA
   while ( in.good() ) { getline(in,ln); if (ln[0] == '!') break; }   

   for (int dim=0; dim<3; ++dim) {
      in >> vel[dim]; getline(in,ln);
   }

}


//------------------------------------------------------------------------------------------------


void frame_name(int mynod, int current_frame_number, char frame[])
{
   int k;
   char proc_id[] = "00", frame_id[] = "0000";
   sprintf(proc_id,"%02d",mynod);
   sprintf(frame_id,"%04d",current_frame_number);
   for (k=0; k<2; ++k) frame[23+k] = proc_id[k];
   for (k=0; k<4; ++k) frame[26+k] = frame_id[k];
}


//------------------------------------------------------------------------------------------------
//transport_output/t_out_00_0000.plt
//transport_output/t_mc0000_00_0000.plt

void realization_frame_name(int mynod, int StochCount, int current_frame_number, char frame[])
{
   int k;
   char mc_id[] = "0000", proc_id[] = "00", frame_id[] = "0000";
   sprintf(mc_id,"%04d",StochCount);
   sprintf(proc_id,"%02d",mynod);
   sprintf(frame_id,"%04d",current_frame_number);
   for (k=0; k<4; ++k) frame[21+k] = mc_id[k];
   for (k=0; k<2; ++k) frame[26+k] = proc_id[k];
   for (k=0; k<4; ++k) frame[29+k] = frame_id[k];
}


//------------------------------------------------------------------------------------------------


void print_output( int& model, eqn *Ptr, double c[][8], fstream& out, int& frame, double& t, int& testnumber, int& mynod)
{
   int i1, i2, i3, p, k;
   int x1index, x2index, x3index;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1*n2*n3;

   double Cerr, Uz;

   int Nd1, Nd2, Nd3, Nd4, Nd5, Nd6, Nd7, Nd8;

   int index_map[8] = {1, 2, 3, 0, 5, 6, 7, 4}; // See below
/*

 transport() uses the following local ordering:

  Y                      3______2
                        /|      /
  |                    / |     /|
  |                   /7_|___6/ |
  |                   |  |    | |
  |______ X           |  0____|_1
  /                   | /     | /
 /                    |/      |/
Z                     /4_____5/



 Tecplot uses the following local ordering:

  Y                      2______1
                        /|      /
  |                    / |     /|
  |                   /6_|___5/ |
  |                   |  |    | |
  |______ X           |  3____|_0
  /                   | /     | /
 /                    |/      |/
Z                     /7_____4/

*/

   out.setf(ios::fixed);
   out.precision(4);

   char proc_id[] = "00";
   sprintf(proc_id,"%02d",mynod);

   out << "VARIABLES = X, Y, Z, C, C_ERROR, U, V, W" << endl;
   out << "ZONE T=\"Proc " << proc_id << "\", ZONETYPE=FEBRICK, DATAPACKING=POINT,   N=" << 8*n << ",   E=" << n << 
          ", STRANDID=" << frame << ", SOLUTIONTIME=" << t << endl;
   out << endl;

   out.unsetf(ios::fixed);
   out.setf(ios::scientific);
   out.precision(15);

   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               x1index = ((p==3)||(p==2)||(p==7)||(p==6)) ? i1 : (i1+1);
               x2index = ((p==3)||(p==0)||(p==7)||(p==4)) ? i2 : (i2+1);
               x3index = ((p==3)||(p==0)||(p==2)||(p==1)) ? i3 : (i3+1);

               out << Ptr->x1[x1index] << "   " << Ptr->x2[x2index] << "   " << Ptr->x3[x3index] << "   ";

               out << c[k][index_map[p]] << "   ";

               Cerr = ( c[k][index_map[p]] -
			   c_true(model, Ptr->x1[x1index], Ptr->x2[x2index], Ptr->x3[x3index], t, testnumber) );
               out << Cerr << "   ";

               Uz = (testnumber==9) ? 0.0 : ( Ptr->u3[k+n*index_map[p]] );            
               out << Ptr->u1[k+n*index_map[p]] << "   " << Ptr->u2[k+n*index_map[p]] << "   " << Uz << endl;
            }
         }
      }
   }

   // connectivity list
   Nd1 = 1; Nd2 = 2; Nd3 = 3; Nd4 = 4; Nd5 = 5; Nd6 = 6; Nd7 = 7; Nd8 = 8;
   for (k=0; k<n; ++k) {
      out << Nd1 << " " << Nd2 << " " << Nd3 << " " << Nd4 << " " << Nd5 << " " << Nd6 << " " << Nd7 << " " << Nd8 << endl;
      Nd1+=8; Nd2+=8; Nd3+=8; Nd4+=8; Nd5+=8; Nd6+=8; Nd7+=8; Nd8+=8;
   }

}


//------------------------------------------------------------------------------------------------


void print_output( eqn *Ptr, double qx[][8], double qy[][8], double qz[][8], fstream& out )
{
   int i1, i2, i3, p, k;
   int x1index, x2index, x3index;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1*n2*n3;

   int Nd1, Nd2, Nd3, Nd4, Nd5, Nd6, Nd7, Nd8;

   int index_map[8] = {1, 2, 3, 0, 5, 6, 7, 4}; // See below
/*

 transport() uses the following local ordering:

  Y                      3______2
                        /|      /
  |                    / |     /|
  |                   /7_|___6/ |
  |                   |  |    | |
  |______ X           |  0____|_1
  /                   | /     | /
 /                    |/      |/
Z                     /4_____5/



 Tecplot uses the following local ordering:

  Y                      2______1
                        /|      /
  |                    / |     /|
  |                   /6_|___5/ |
  |                   |  |    | |
  |______ X           |  3____|_0
  /                   | /     | /
 /                    |/      |/
Z                     /7_____4/

*/

   out.setf(ios::scientific);
   out.precision(15);
   out << "VARIABLES = X, Y, Z, Qx, Qy, Qz, U, V, W" << endl;
   out << "ZONE   ZONETYPE=FEBRICK, DATAPACKING=POINT,   N=" << 8*n << ",   E=" << n << endl;
   out << endl;

   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               x1index = ((p==3)||(p==2)||(p==7)||(p==6)) ? i1 : (i1+1);
               x2index = ((p==3)||(p==0)||(p==7)||(p==4)) ? i2 : (i2+1);
               x3index = ((p==3)||(p==0)||(p==2)||(p==1)) ? i3 : (i3+1);

               out << Ptr->x1[x1index] << "   " << Ptr->x2[x2index] << "   " << Ptr->x3[x3index] << "   ";
               out << qx[k][index_map[p]] << "   " << qy[k][index_map[p]] << "   " << qz[k][index_map[p]] << "   ";
               out << Ptr->u1[k+n*index_map[p]] << "   " << Ptr->u2[k+n*index_map[p]] << "   " << Ptr->u3[k+n*index_map[p]] << endl;
            }
         }
      }
   }

   // connectivity list
   Nd1 = 1; Nd2 = 2; Nd3 = 3; Nd4 = 4; Nd5 = 5; Nd6 = 6; Nd7 = 7; Nd8 = 8;
   for (k=0; k<n; ++k) {
      out << Nd1 << " " << Nd2 << " " << Nd3 << " " << Nd4 << " " << Nd5 << " " << Nd6 << " " << Nd7 << " " << Nd8 << endl;
      Nd1+=8; Nd2+=8; Nd3+=8; Nd4+=8; Nd5+=8; Nd6+=8; Nd7+=8; Nd8+=8;
   }

}

//------------------------------------------------------------------------------------------------

void write_binary( eqn *Ptr, double cexp[][8], double cvar[][8], fstream& out)
{
   int i1, i2, i3, p, k;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int index_map[8] = {1, 2, 3, 0, 5, 6, 7, 4}; // See above

   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               out.write((char *)(&cexp[k][index_map[p]]), sizeof(double));
            }
         }
      }
   }
   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               out.write((char *)(&cvar[k][index_map[p]]), sizeof(double));
            }
         }
      }
   }
}



//------------------------------------------------------------------------------------------------
// miz17, bag8  6/2/10
// read_add_stochastic performs the following tasks:
//   read in current cexp and cvar from fstream& in
//   compute cexp = cexp + cweight*c
//   compute cvar = cvar + cweight*(c**2)
//   then return new cexp and cvar


// ***************************************
// miz17, bag8 - 6/29/10
// This is now binary!!!!!!!
// ***************************************

void read_add_stochastic( eqn *Ptr, double c[][8], double cexp[][8], double cvar[][8],
                          double &cweight, fstream& in, bool first_realization, int& mynod)
{
   int i1, i2, i3, p, k;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int index_map[8] = {1, 2, 3, 0, 5, 6, 7, 4}; // See above

   // read in cexp and cvar, and compute new values
   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               if (! first_realization) {
                  in.read((char *)(&cexp[k][index_map[p]]), sizeof(double));
                  cexp[k][index_map[p]] = cexp[k][index_map[p]] +
                     cweight * c[k][index_map[p]];
               }
               else {
                  cexp[k][index_map[p]] = cweight * c[k][index_map[p]];
               }
            }
         }
      }
   }
   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               if (! first_realization) {
                  in.read((char *)(&cvar[k][index_map[p]]), sizeof(double));
                  cvar[k][index_map[p]] = cvar[k][index_map[p]] +
                     cweight * pow(c[k][index_map[p]],2) ;
               }
               else {
                  cvar[k][index_map[p]] = cweight * pow(c[k][index_map[p]],2);
               }
            }
         }
      }
   }
}


//------------------------------------------------------------------------------------------------


void subtract_square_exp( eqn *Ptr, double cexp[][8], double cvar[][8])
{

   int i1, i2, i3, p, k;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1*n2*n3;
   int index_map[8] = {1, 2, 3, 0, 5, 6, 7, 4}; // See above

   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               cvar[k][index_map[p]] = cvar[k][index_map[p]] - pow(cexp[k][index_map[p]],2);
            }
         }
      }
   }
}


//------------------------------------------------------------------------------------------------


void print_stochastic( eqn *Ptr, double cexp[][8], double cvar[][8], fstream& out, int& frame, double& t, int& mynod)
{

  int i1, i2, i3, p, k;
   int x1index, x2index, x3index;
   int n1 = Ptr->n1, n2 = Ptr->n2, n3 = Ptr->n3;
   int n = n1*n2*n3;
 
   double Cerr, Uz;

   int Nd1, Nd2, Nd3, Nd4, Nd5, Nd6, Nd7, Nd8;

   int index_map[8] = {1, 2, 3, 0, 5, 6, 7, 4}; // See above

   out.setf(ios::fixed);
   out.precision(4);

   char proc_id[] = "00";
   sprintf(proc_id,"%02d",mynod);

   out << "VARIABLES = X, Y, Z, CEXP, CVAR" << endl;
   out << "ZONE T=\"Proc " << proc_id << "\", ZONETYPE=FEBRICK, DATAPACKING=POINT,   N=" << 8*n << ",   E=" << n << 
          ", STRANDID=" << frame << ", SOLUTIONTIME=" << t << endl;
   out << endl;

   out.unsetf(ios::fixed);
   out.setf(ios::scientific);
   out.precision(15);

   for (i3=0; i3<n3; ++i3) {
      for (i2=0; i2<n2; ++i2) {
         for (i1=0; i1<n1; ++i1) {
            k = i1 + n1*i2 + n1*n2*i3;
            for (p=0; p<8; ++p) {
               x1index = ((p==3)||(p==2)||(p==7)||(p==6)) ? i1 : (i1+1);
               x2index = ((p==3)||(p==0)||(p==7)||(p==4)) ? i2 : (i2+1);
               x3index = ((p==3)||(p==0)||(p==2)||(p==1)) ? i3 : (i3+1);

               out << Ptr->x1[x1index] << "   " << Ptr->x2[x2index] << "   " << Ptr->x3[x3index] << "   ";

               out << cexp[k][index_map[p]] << "   ";
               
               out << cvar[k][index_map[p]] << endl;
            }
         }
      }
   }

   // connectivity list
   Nd1 = 1; Nd2 = 2; Nd3 = 3; Nd4 = 4; Nd5 = 5; Nd6 = 6; Nd7 = 7; Nd8 = 8;
   for (k=0; k<n; ++k) {
      out << Nd1 << " " << Nd2 << " " << Nd3 << " " << Nd4 << " " << Nd5 << " " << Nd6 << " " << Nd7 << " " << Nd8 << endl;
      Nd1+=8; Nd2+=8; Nd3+=8; Nd4+=8; Nd5+=8; Nd6+=8; Nd7+=8; Nd8+=8;
   }

}


//------------------------------------------------------------------------------------------------
