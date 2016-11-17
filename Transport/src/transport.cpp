#include "transport.h"

void transport_(int& mynod, int& model, int *mynbr, int &n1, int &n2, int &n3,
		double *u1, double *u2, double *u3,
		double *x1, double *x2, double *x3,
		int &itest, int &StochFlag, double &cweight,
        int *RealizationIdx, int *NumOfRealizations, int *totalloops)
// mynod             :   processor ID
// model             :   subdomain model, 1=Darcy or 2=Stokes
// mynbr             :   IDs of the neighbour processors

// n1, n2, n2        :   number of divisions in the local subdomain
//                       along the corresponding coordinate axis

// u1, u2, u3        :   components of the velocity
// x1, x2, x3        :   arrays of length (n? + 1) to hold the coordinates of the nodes

// itest             :   test number

// StochFlag         :   if StochFlag=1 KL expansion for the permeability is used for the porous medium
//                       (many realizations)
//
//                       if StochFlag=0 deterministic permeability for the porous medium or pure Stokes flow
//                       (one realization)

// cweight           :   used in the stochastic loop to compute the expectation and the variance

// RealizationIdx    :   array containing the indices of stochastic realizations to plot 
//                       (if KL expansion is used for the permeability in the porous medium)

// NumOfRealizations :   length of RealizationIdx

{

   static int StochCount = 0; // how many times the function is called

   //---- THE INPUT DATA IS READ BY PROCESSOR 0 AND BROADCASTED TO THE REST ----
   int tflag, velflag, testnumber, IntegrationMethod, limflag, ldgerr, dflag, dtype;
   double bdry[3][2], d_molecular, d_longitudinal, d_transverse, porosity, d_visc, viscosity;
   double D[9], T, dt, alpha;               // diffusion, time, step, limiter parameter
   int skip;                                // number of steps to skip before making a frame
   
   int *mpi_i_buffer    = new int    [10];   // auxiliary variables, which are used
   double *mpi_d_buffer = new double [23];  // for message broadcasting
  
   if (mynod==0) {
      ifstream in("input"); // open the input file
      read_file(in, bdry, tflag, velflag, testnumber, limflag, alpha,
                ldgerr, dflag, dtype, D, d_molecular, d_longitudinal, d_transverse, 
                porosity, d_visc, IntegrationMethod, T, dt, skip);
      in.close();  // close the input file


      mpi_i_buffer[0] = tflag;
      mpi_i_buffer[1] = velflag;
      mpi_i_buffer[2] = testnumber;
      mpi_i_buffer[3] = itest;
      mpi_i_buffer[4] = limflag;
      mpi_i_buffer[5] = ldgerr;
      mpi_i_buffer[6] = dflag;
      mpi_i_buffer[7] = dtype;
      mpi_i_buffer[8] = skip;
      mpi_i_buffer[9] = IntegrationMethod;

      mpi_d_buffer[0] = bdry[0][0]; mpi_d_buffer[1] = bdry[0][1];
      mpi_d_buffer[2] = bdry[1][0]; mpi_d_buffer[3] = bdry[1][1];
      mpi_d_buffer[4] = bdry[2][0]; mpi_d_buffer[5] = bdry[2][1];
      mpi_d_buffer[6] = alpha; mpi_d_buffer[7] = T; mpi_d_buffer[8] = dt;
      for (int i=0; i<9; ++i) mpi_d_buffer[i+9] = D[i];
      mpi_d_buffer[18] = d_molecular;
      mpi_d_buffer[19] = d_longitudinal;
      mpi_d_buffer[20] = d_transverse;
      mpi_d_buffer[21] = porosity;
      mpi_d_buffer[22] = d_visc;
   }
   gibcastX(4000, mpi_i_buffer, 10, 0);
   gdbcastX(5000, mpi_d_buffer, 23, 0);
   if (mynod>0) {
      tflag      = mpi_i_buffer[0];
      velflag    = mpi_i_buffer[1];
      testnumber = mpi_i_buffer[2];
      itest      = mpi_i_buffer[3];
      limflag    = mpi_i_buffer[4];
      ldgerr     = mpi_i_buffer[5];
      dflag      = mpi_i_buffer[6];
      dtype      = mpi_i_buffer[7];
      skip       = mpi_i_buffer[8];
      IntegrationMethod = mpi_i_buffer[9];

      bdry[0][0] = mpi_d_buffer[0]; bdry[0][1] = mpi_d_buffer[1];
      bdry[1][0] = mpi_d_buffer[2]; bdry[1][1] = mpi_d_buffer[3];
      bdry[2][0] = mpi_d_buffer[4]; bdry[2][1] = mpi_d_buffer[5];
      alpha = mpi_d_buffer[6]; T = mpi_d_buffer[7]; dt = mpi_d_buffer[8];
      for (int i=0; i<9; ++i) D[i] = mpi_d_buffer[i+9];
      d_molecular    = mpi_d_buffer[18];
      d_longitudinal = mpi_d_buffer[19];
      d_transverse   = mpi_d_buffer[20];
      porosity       = mpi_d_buffer[21];
      d_visc         = mpi_d_buffer[22];
   }
   delete mpi_i_buffer;
   delete mpi_d_buffer;
   //---------------------------------------------------------------------------   

   if (mynod==0) system("mkdir -p transport_output");
   gsyncX();

   //-------------- create a unique filename for each processor ----------------
   //               ( used for debugging)


   //  the filename string is manipulated using C++ strings 
   //string name = "transport_output/dbg", proc_id = "00";
   //sprintf(&proc_id[0],"%02d",mynod);
   //name = name + "_" + proc_id + ".plt";
   //const char *filename = name.c_str();

   //   the filename string is manipulated in a C like manner
   char filename[] = "transport_output/dbg_00.plt", proc_id[] = "00";
   sprintf(proc_id,"%02d",mynod);
   filename[21] = proc_id[0]; filename[22] = proc_id[1];

   //---------------------------------------------------------------------------

   //------------------- open a debug file for each processor ------------------
   fstream out;
   out.open( filename, ios::out );
   out.setf(ios::fixed);
   //---------------------------------------------------------------------------


   eqn eqn_data(n1, n2, n3, bdry, x1, x2, x3, u1, u2, u3), *Ptr = &eqn_data;

   double (*c)[8] = new double [n1*n2*n3][8];  // array to store the initial values of the concentration
                                               // Note: the local ordering is explained in utilities.h

   // miz17, bag8 - declare concentration expectation and variance arrays
   double (*cexp)[8], (*cvar)[8];
   if (StochFlag==1) {
      cexp = new double [n1*n2*n3][8];
      cvar = new double [n1*n2*n3][8];
   }
   
   /********** USE IT TO READ DEBUG DATA FROM THE INPUT FILE *********
   // currently reads just the components of the velocity
   double vel[3];
   ifstream in_dbg("input");
   read_file_dbg(in_dbg, vel);
   in_dbg.close();
   set_velocity(Ptr, vel); // set the velocity to be a constant vector 
   *******************************************************************/

   if (  ( ((testnumber>=1)&&(testnumber<=12))||(testnumber==14) )&&(dtype!=1)  ) {
      cout << "\ntransport test " << testnumber << " needs dtype=1\n" << endl;
      exit(1);
   }
   if(  ( (testnumber==7)||(testnumber==8)||(testnumber==13)||
          (testnumber==15)||(testnumber==16) )&&(ldgerr==1)  ) {
      cout << "\ntransport test " << testnumber << " needs ldgerr=0\n" << endl;
      exit(1);
   }
   if(  ( (testnumber==2)||(testnumber==4)||(testnumber==6)||(testnumber==8) )&&( (dflag!=1) )  ) {
      cout << "\ntransport test " << testnumber << " needs dflag=1\n" << endl;
      exit(1);
   }
   if(  ( (testnumber==1)||(testnumber==3)||(testnumber==5)||(testnumber==7) )&&(dflag!=0)  ) {
      cout << "\ntransport test " << testnumber << " needs dflag=0\n" << endl;
      exit(1);
   }

   if (  ( ((testnumber>=9)&&(testnumber<=12))||(testnumber==14) )&&(dflag==0)  ) {
      cout << "\ntransport test " << testnumber << " needs dflag=1\n";
      cout << "To run this test without diffusion\n";
      cout << "assign zeros to the D entries in the input file.\n" << endl;
      exit(1); 
   }

   if( (testnumber==10)&&(itest!=-10) ) {
      cout << "\ntransport test " << testnumber << " needs the SD test to be -10\n" << endl;
      exit(1);
   }
   if( (testnumber==11)&&(itest!=-11) ) {
      cout << "\ntransport test " << testnumber << " needs the SD test to be -11\n" << endl;
      exit(1);
   }
   if( (testnumber==12)&&(itest!=-12) ) {
      cout << "\ntransport test " << testnumber << " needs the SD test to be -12\n" << endl;
      exit(1);
   }
   if( (testnumber==14)&&(itest!=-100)&&(itest!=2150) ) {
      cout << "\ntransport test " << testnumber << " needs the SD test to be -100 or 2150\n" << endl;
      cout << "testnumber = " << testnumber << endl;
      cout << "itest = " << itest << endl;
      exit(1);
   }

   //---------------------------------- SET THE DIFFUSION MATRIX ------------------------------
   if (dtype==1) {
      for (int i=0; i<3; ++i) {
	 for (int j=0; j<3; ++j) {
	    Ptr->D[i][j] = D[i+3*j];
	 }
      }
   }

   else if (dtype==2) {

      if (model==2) {
         int param = 6; viscosity = sdparameters_(&param); // get the viscosity
         Ptr->D_fluid = viscosity * d_visc;
      }

      if (model==1) {
         double unorm, uavr[3];
         int nel = n1*n2*n3;
	 Ptr->D_porous = new double [nel][3][3];
	 for (int i3=0; i3<n3; ++i3) {
	    for (int i2=0; i2<n2; ++i2) {
	       for (int i1=0; i1<n1; ++i1) {
		  
		  int k = i1 + n1*i2 + n1*n2*i3; // element
	       
		  // use the cell averaged velocity
		  uavr[0] = 0.0; uavr[1] = 0.0; uavr[2] = 0.0;
		  for (int p=0; p<8; ++p) {
		     uavr[0] +=  u1[k+p*nel]; uavr[1] +=  u2[k+p*nel]; uavr[2] +=  u3[k+p*nel]; 
		  }
		  uavr[0] = uavr[0]/8.0; uavr[1] = uavr[1]/8.0; uavr[2] = uavr[2]/8.0;
		  unorm = sqrt(uavr[0]*uavr[0] + uavr[1]*uavr[1] + uavr[2]*uavr[2]);
		  for (int i=0; i<3; ++i) {
		     for (int j=0; j<3; ++j) {
			double t1 = (i==j) ? ( porosity*d_molecular ) : 0.0;
			double t2 = d_longitudinal * uavr[i]*uavr[j] / unorm;
			double t3 = (i==j) ? ( unorm*d_transverse ) : 0.0;
			double t4 = d_transverse * uavr[i]*uavr[j] /unorm;                        
			Ptr->D_porous[k][i][j] = t1 + t2 + t3 - t4;
		     }
		  }
		  
	       }
	    }
	 }
	 
      }

   }
   //------------------------------------------------------------------------------------------


   //--------------------------------- SET THE POROSITY ---------------------------------------
   //  for model=2 (Stokes) the default value for the porosity is 1.0

   //  leave the old test cases unaffected
   if (  ( (testnumber==9)||(testnumber>14) ) && (model==1)  ) Ptr->porosity  = porosity; 
   //------------------------------------------------------------------------------------------


   //----- use Cholesky factorization for the mass matrix -----
   double *temporary = new double [8];
   for (int i=0; i<64; ++i) Ptr->M[i] = Ptr->M[i] * Ptr->porosity; 
   char uplo='U';
   int dim=8, nrhs=1, info;
   dposv_(uplo, dim, nrhs, Ptr->M, dim, temporary, dim, info);
   delete [] temporary;
   //----------------------------------------------------------

         
   if (mynod==0) {

      ofstream param_out("transport_output/parameters");
      param_out.setf(ios::fixed);
      param_out << "\n\nPARAMETERS FOR THE TRANSPORT CODE\n" << endl;
      param_out << "tflag      = " << tflag << endl;
      param_out << "velflag    = " << velflag << endl;
      param_out << "testnumber = " << testnumber << endl;
      param_out << "\n d_molecular    = " << d_molecular << endl;
      param_out << " d_longitudinal = " << d_longitudinal << endl;
      param_out << " d_transverse   = " << d_transverse << endl;
      param_out << " porosity       = " << Ptr->porosity << endl;
      param_out << " d_visc         = " << d_visc << endl;
      param_out << " viscosity      = " << viscosity << endl;
      param_out << " limflag        = " << limflag << endl;
      param_out << " ldgerr         = " << ldgerr << endl;
      param_out << " dflag          = " << dflag << endl;
      param_out << " dtype          = " << dtype << endl;
      param_out << " integration method   = " << IntegrationMethod << endl;
      param_out << " T = " << T << "   dt = " << dt << endl;
      param_out << " frames to be skipped = " << skip << endl;
      param_out << "\ndomain boundary:" << endl;
      param_out << bdry[0][0] << "   " << bdry[0][1] << endl;
      param_out << bdry[1][0] << "   " << bdry[1][1] << endl;
      param_out << bdry[2][0] << "   " << bdry[2][1] << endl << endl;
      if ( (dflag!=0)&&(dtype==1) ) for (int i=0; i<3; ++i) param_out << "D[" << i << "][0] = " 
                                                                      << Ptr->D[i][0]
								      << "   D[" << i << "][1] = "
								      << Ptr->D[i][1]
								      << "   D[" << i << "][2] = "
								      << Ptr->D[i][2] << endl;
      if ( (0<testnumber)&&(testnumber<8) ) param_out << "Ux = " << Ux << endl << endl;
      param_out << endl;
      param_out.close();

      /***************************************************************
      cout << "\n";
      //for (int i=0; i<3; ++i) cout << "vel[" << i << "] = " << vel[i] << endl;
      //cout << "\n";
      cout << "\n\n" << endl;
      ****************************************************************/

   }

   /*
   function w[4];

   w[0] = w0;
   w[1] = w1;
   w[2] = w2;
   w[3] = w3;
   for (int i=0; i<4; ++i) {
      cout << quad(2, -1.0, w[i], w[0]) << " " << quad(2, -1.0, w[i], w[1]) << " "
           << quad(2, -1.0, w[i], w[2]) << " " << quad(2, -1.0, w[i], w[3]) << endl;
   }
   cout << "----------" << endl;

   w[0] = w0;
   w[1] = w3;
   w[2] = w7;
   w[3] = w4;
   for (int i=0; i<4; ++i) {
      cout << quad(0, -1.0, w[i], w[0]) << " " << quad(0, -1.0, w[i], w[1]) << " "
           << quad(0, -1.0, w[i], w[2]) << " " << quad(0, -1.0, w[i], w[3]) << endl;
   }
   cout << "----------" << endl;
   w[0] = w0;
   w[1] = w4;
   w[2] = w5;
   w[3] = w1;
   for (int i=0; i<4; ++i) {
      cout << quad(1, -1.0, w[i], w[0]) << " " << quad(1, -1.0, w[i], w[1]) << " "
           << quad(1, -1.0, w[i], w[2]) << " " << quad(1, -1.0, w[i], w[3]) << endl;
   }
   */

    
   // initial condition
   apply_IC_discts(model, n1, x1, n2, x2, n3, x3, c, testnumber);
   //apply_IC(model, n1, x1, n2, x2, n3, x3, c, testnumber);

   // explicit time stepping
   switch (IntegrationMethod) {
   case 1 :
      fwd_euler( mynod, model, mynbr, Ptr, limflag, alpha, dflag, dtype, c, cexp, cvar, 
                 dt, T, skip, out, testnumber, ldgerr, StochFlag, StochCount, cweight,
                 RealizationIdx, NumOfRealizations, totalloops );
      break;
   case 2 :
      rk2( mynod, model, mynbr, Ptr, limflag, alpha, dflag, dtype, c, cexp, cvar, 
           dt, T, skip, out, testnumber, ldgerr, StochFlag, StochCount, cweight,
           RealizationIdx, NumOfRealizations, totalloops );
      break;
   case 4 :
      rk4( mynod, model, mynbr, Ptr, limflag, alpha, dflag, dtype, c, cexp, cvar,
           dt, T, skip, out, testnumber, ldgerr, StochFlag, StochCount, cweight,
           RealizationIdx, NumOfRealizations, totalloops );
      break;
   default :
      cout << "\nThe integration method " << IntegrationMethod << " is not implemented!\n" << endl;
      exit(1);
   }

   ++StochCount;

   out.close(); // close the debug otput file

   if ( (dtype==2)&&(model==1) ) delete [] (Ptr->D_porous);

   delete [] c;
   if (StochFlag==1) {
      delete [] cexp;
      delete [] cvar;
   }

}
