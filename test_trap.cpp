/************************************************************************
 * File: test_trap.cpp                                              *
 * Programmers: Cory Mikida                                             *     
 *----------------------------------------------------------------------*
 * Testing hand-coded trapezoidal integration of chemistry with PyJac generated source
 * terms, with and without Cantera involvement for dependent variable
 * update.
 * Example problem: homogeneous Cantera mechanism (9-species San-Diego) *
 ************************************************************************/

#include "jacobian.hpp"
#include "species_rates.hpp"
#include "memcpy_2d.hpp"
#include <stdio.h>
#include <iostream>     // std::cout, std::fixed
#include <iomanip>
#include <fstream>

/* Problem Constants */

#define NSP   9             /* number of species */
#define NEQ   10            /* number of equations  */
#define RTOL  1e-6         /* scalar relative tolerance            */
#define ATOL 1e-12         /* absolute tolerance component */
#define T0    0.0          /* initial time           */
#define T1    2.5e-8          /* first output time      */
#define TADD  2.5e-8         /* output time factor     */
#define NOUT  3800           /* number of output times */
#define USE_CANTERA 1


/***************************** Main Program ******************************/

using namespace std;

extern "C" void dgesv_(int *, int *, double *, int *, int *, double *, int *, int * );
extern "C" double dnrm2_(int *, double *, int *);

main()
{
  double reltol, abstol, tout;
  int iout, flag;
  double y[NSP+1];
  double mw[NSP];
  double corr_norm;

  // For V2: need number of moles of each species as input.
  y[0] = 1000.0;            /* PyJacV2 convention: temperature goes at the beginning. */
  //y[1] = 101325.0;            /* PyJacV2 convention: pressure goes at the beginning. */
  y[1] = 101323.2510934631;            /* PyJacV2 convention: pressure goes at the beginning. */
  y[2] = 0.0;               /* Leave out last species as well */
  y[3] = 0.0;
  y[4] = (0.1949817389200555/0.2072648773462248)/31.998;
  y[5] = (0.01228313842616932/0.2072648773462248)/2.016;
  y[6] = 0.0;
  y[7] = 0.0;
  y[8] = 0.0;
  y[9] = 0.0;

  mw[0] = 15.9994;
  mw[1] = 1.00794;
  mw[2] = 15.9994*2;
  mw[3] = 1.00794*2;
  mw[4] = 15.9994 + 1.00794;
  mw[5] = 15.9994 + 2*1.00794;
  mw[6] = 2*15.9994 + 1.00794;
  mw[7] = 2*15.9994 + 2*1.00794;
  mw[8] = 2*14.00674;
 
  double R = 8314.4621; /* gas constant */  
  double pres = 101325; /* Pressure in pa */
  double rho = 0.2072648773462248;  /* density */
  double vol = 1.0 / 0.2072648773462248;  /* specific volume */

  /* Initialize jacobian to pass to generated code */
  double jac[(NSP+1)*(NSP+1)];
  double jac_trans[(NSP+1)*(NSP+1)];
  /* Initialize source to pass to generated code */
  double dy[NSP+1];
  double dy_old[NSP+1];

  /* Work arrays */
  double* rwk_dy = (double*)malloc(245 * sizeof(double));
  memset(rwk_dy, 0, 245 * sizeof(double));
  double* rwk_jac = (double*)malloc(363 * sizeof(double));
  memset(rwk_jac, 0, 363 * sizeof(double));

  reltol = RTOL;               /* Set the scalar relative tolerance */
  abstol = ATOL;               /* Set the scalar absolute tolerance */

  /* For Lapack */
  int ipiv[NSP+1], info;
  int nrhs = 1;
  double yguess[NSP+1];
  double yold[NSP+1];
  double corr[NSP+1];
  double corr_weights[NSP+1];
  double corr_weighted[NSP+1];
  int nsp_l = NSP+1;

  /* 1D identity matrix */
  double ident[(NSP+1)*(NSP+1)];
  for (int i=0; i < (NSP+1)*(NSP+1); i++) {
    if (i % (NSP + 2) == 0) {
      ident[i] = 1;
    } else {
      ident[i] = 0;
    } 
  }

  for (int j=0;j <= NSP; j++) {
    yold[j] = y[j];
  }
  /* Newton iteration using Jacobian, printing results*/
  printf(" \nHomogeneous Cantera problem\n\n");

  // Remove old copies of output files
  for (int i=2; i < 10; i++) {
    if (remove(("Output/test_trap_" + std::to_string(i-2) + ".txt").c_str()) != 0) {
      perror ("error deleting file");
    }
  }
  if (remove("Output/test_trap_temperature.txt") != 0) {
    perror ("error deleting file");
  }

  /* Timestepping loop */
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += TADD) {
    /* Simple Newton loop. */
    corr_norm = 1.0;
    for (int j=0;j <= NSP; j++) {
      yguess[j] = y[j];
    }
    double tin = tout - TADD;
    species_rates (&tin, &vol, y, dy_old, rwk_dy);  /* source term */
    while (abs(corr_norm) >= reltol) {
    //while (abs(corr_norm) >= 1e-9) {
      species_rates (&tout, &vol, yguess, dy, rwk_dy);  /* source term */
      jacobian (&tout, &vol, yguess, jac, rwk_jac); /* Jacobian evaluation */
      for (int j=0;j <= NSP; j++) {
        dy[j] = yguess[j] - yold[j] - 0.5 * TADD * (dy_old[j] + dy[j]); // Trap
      }

      // Transpose the Jacobian
      // PYJAC outputs Fortran-ordering
      for (int i = 0; i < (NSP+1); ++i )
      {
         for (int j = 0; j < (NSP+1); ++j )
         {
            // Index in the original matrix.
            int index1 = i*(NSP+1)+j;

            // Index in the transpose matrix.
            int index2 = j*(NSP+1)+i;

            jac_trans[index2] = jac[index1];
         }
      }

      for (int i=0; i<(NSP+1)*(NSP+1); i++) {
          jac[i] = jac_trans[i];
      }
      /* Subtract from identity to get jac for Newton */
      for (int j=0;j < (NSP+1)*(NSP+1); j++) {
        jac[j] = ident[j] - TADD * 0.5 * jac[j]; // Trap
        jac[j] = -jac[j];
      }
      /* Call LAPACK to invert */
      dgesv_(&nsp_l, &nrhs, jac, &nsp_l, ipiv, dy, &nsp_l, &info);
      for (int j=0;j <= NSP; j++) {
        corr[j] = dy[j];
        corr_weights[j] = 1.0 / (reltol * abs(yguess[j]) + abstol);
	corr_weighted[j] = corr[j]*corr_weights[j];
      }
      // FIXME for now zero out the pressure
      corr[1] = 0.0;
      corr_weighted[1] = 0.0;
      //corr_norm = dnrm2_(&nsp_l, corr, &nrhs);
      corr_norm = dnrm2_(&nsp_l, corr_weighted, &nrhs);
      for (int j=0;j <= NSP; j++) {
        yguess[j] = yguess[j] + corr[j];
      }
      // FIXME: try self pressure calc.
      double rho_test = 0.0;
      for (int j=2;j < NSP; j++) {
        rho_test += yguess[j]*rho;
      }
      pres = rho_test*yguess[0]*R;
      yguess[1] = pres;
    }
      
    for (int j=0;j <= NSP; j++) {
      y[j] = yguess[j];
      yold[j] = y[j];
    }

    printf("At t = %0.4e      y =%14.6e\n",
           tout, y[0]);

    // Print solutions to file
    for (int i=2; i < 10; i++) {
      std::ofstream outfile;
      outfile.open("Output/test_trap_" + std::to_string(i-2) + ".txt", std::ios_base::app);
      outfile << std::setprecision(16) << y[i] << std::endl;
    }
    std::ofstream outfile;
    outfile.open("Output/test_trap_temperature.txt", std::ios_base::app);
    outfile << std::setprecision(16) << y[0] << std::endl;
  }

  return(0);
}
