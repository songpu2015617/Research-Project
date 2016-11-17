#include "LDG.h"

extern "C" {

  void transport_(int &mynod, int &model, int *mynbr, int &n1, int &n2, int &n3,
		  double *u1, double *u2, double *u3,
		  double *x1, double *x2, double *x3,
		  int &itest, int &StochFlag, double &cweight,
		  int *RealizationIdx, int *NumOfRealizations, int *totalloops);

// mynod        :   processor ID
// model        :   subdomain model, 1=Darcy or 2=Stokes
// mynbr        :   IDs of the neighbour processors

// n1, n2, n2   :   number of divisions in the local subdomain
//                  along the corresponding coordinate axis

// u1, u2, u3   :   components of the velocity
// x1, x2, x3   :   arrays of length (n? + 1) to hold the coordinates of the nodes

// itest        :   test number

// StochFlag    :   if StochFlag=1 KL expansion for the permeability is used for the porous medium
//                  (many realizations)
//
//                  if StochFlag=0 deterministic permeability for the porous medium or pure Stokes flow
//                  (one realization)

// cweight      :   used in the stochastic loop to compute the expectation and the variance

};



