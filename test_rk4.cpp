/************************************************************************
 * File: test_rk4.cpp                                                   *
 * Programmers: Cory Mikida                                             *     
 *----------------------------------------------------------------------*
 * Testing RK4 integration of chemistry with PyJac generated source
 * terms, with and without Cantera involvement for dependent variable
 * update.
 * Example problem: homogeneous 9-species San-Diego mechanism *
 ************************************************************************/

#include "cantera/IdealGasMix.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include "jacobian.hpp"
#include "species_rates.hpp"
#include "memcpy_2d.hpp"
#include <stdio.h>
#include <iostream>
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
#define USE_CANTERA 0


/***************************** Main Program ******************************/

using namespace std;

main()
{
  double reltol, abstol, tout, tin;
  int iout, flag;
  double y[NSP+1];
  double y1[NSP+1];
  double y2[NSP+1];
  double y3[NSP+1];
  double y4[NSP+1];
  double k1[NSP+1];
  double k2[NSP+1];
  double k3[NSP+1];
  double k4[NSP+1];
  double mw[NSP];
  double massFractions[NSP];

  // PyJac takes number of moles of each species as input.
  y[0] = 1001.0;            /* PyJacV2 convention: temperature goes at the beginning. */
  //y[1] = 101325.0;            /* PyJacV2 convention: pressure goes at the beginning. */
  y[1] = 101426.325;            /* PyJacV2 convention: pressure goes at the beginning. */
  y[2] = 0.0;               /* Leave out last species as well */
  y[3] = 0.0;
  y[4] = (0.1949817389200555/0.2072648773462248)/(2*15.9994);
  y[5] = (0.01228313842616932/0.2072648773462248)/(2*1.00794);
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

  /* Work arrays */
  double* rwk_dy = (double*)malloc(245 * sizeof(double));
  memset(rwk_dy, 0, 245 * sizeof(double));

  /* Newton iteration using Jacobian, printing results*/
  printf(" \nHomogeneous Cantera problem\n\n");

  Cantera::IdealGasMix * GasMixture;
  GasMixture = new Cantera::IdealGasMix("Mechanisms/sanDiego.xml");
  for (int i=2; i < 10; i++) {
    massFractions[i-2] = y[i]*mw[i-2];
  }
  massFractions[8] = 0.0;
  GasMixture->setState_TPY(y[0], y[1], massFractions);
  double int_energy = GasMixture->intEnergy_mass();
  std::cout << "Internal energy from Cantera of IC: " << std::setprecision(16) << int_energy << std::endl;
  double tol = 1e-10;

  // Remove old copies of output files
  for (int i=2; i < 10; i++) {
    if (remove(("Output/test_rk4_" + std::to_string(i-2) + ".txt").c_str()) != 0) {
      perror ("error deleting file");
    }
  }
  if (remove("Output/test_rk4_temperature.txt") != 0) {
    perror ("error deleting file");
  }

  /* Timestepping loop */
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += TADD) {
    /* RK4 timestepper */
    /* STAGE 1 */
    tin = tout - TADD;
    for (int j=0;j <= NSP; j++) {
      y1[j] = y[j];
    }
    if (USE_CANTERA == 1) {
      // Update pressure and temperature using Cantera + constant internal energy assumption.
      for (int i=2; i < 10; i++) {
        massFractions[i-2] = y1[i]*mw[i-2];
      }
      massFractions[8] = 0.0;
      GasMixture->setMassFractions(massFractions);
      GasMixture->setState_UV(int_energy, vol, tol);
      y1[0] = GasMixture->temperature();
      y1[1] = GasMixture->pressure();
    }
    species_rates (&tin, &vol, y1, k1, rwk_dy);  /* source term */
    
    /* STAGE 2 */
    tin = tout - 0.5*TADD;
    for (int j=0;j <= NSP; j++) {
      y2[j] = y1[j] + k1[j]*TADD*0.5;
    }
    if (USE_CANTERA == 1) {
      // Update pressure and temperature using Cantera + constant internal energy assumption.
      for (int i=2; i < 10; i++) {
        massFractions[i-2] = y2[i]*mw[i-2];
      }
      massFractions[8] = 0.0;
      GasMixture->setMassFractions(massFractions);
      GasMixture->setState_UV(int_energy, vol, tol);
      y2[0] = GasMixture->temperature();
      y2[1] = GasMixture->pressure();
    }
    species_rates (&tout, &vol, y2, k2, rwk_dy);  /* source term */

    /* STAGE 3 */
    tin = tout - 0.5*TADD;
    for (int j=0;j <= NSP; j++) {
      y3[j] = y1[j] + k2[j]*TADD*0.5;
    }
    if (USE_CANTERA == 1) {
      // Update pressure and temperature using Cantera + constant internal energy assumption.
      for (int i=2; i < 10; i++) {
        massFractions[i-2] = y3[i]*mw[i-2];
      }
      massFractions[8] = 0.0;
      GasMixture->setMassFractions(massFractions);
      GasMixture->setState_UV(int_energy, vol, tol);
      y3[0] = GasMixture->temperature();
      y3[1] = GasMixture->pressure();
    }
    species_rates (&tout, &vol, y3, k3, rwk_dy);  /* source term */

    /* STAGE 4 */
    tin = tout;
    for (int j=0;j <= NSP; j++) {
      y4[j] = y1[j] + k3[j]*TADD;
    }
    if (USE_CANTERA == 1) {
      // Update pressure and temperature using Cantera + constant internal energy assumption.
      for (int i=2; i < 10; i++) {
        massFractions[i-2] = y4[i]*mw[i-2];
      }
      massFractions[8] = 0.0;
      GasMixture->setMassFractions(massFractions);
      GasMixture->setState_UV(int_energy, vol, tol);
      y4[0] = GasMixture->temperature();
      y4[1] = GasMixture->pressure();
    }
    species_rates (&tout, &vol, y4, k4, rwk_dy);  /* source term */
      
    /* Add 'em up */
    for (int j=0;j <= NSP; j++) {
      y[j] = y[j] + TADD*((1.0/6.0)*k1[j] + (1.0/3.0)*k2[j] + (1.0/3.0)*k3[j] + (1.0/6.0)*k4[j]);
    }
    if (USE_CANTERA == 1) {
      // Update pressure and temperature using Cantera + constant internal energy assumption.
      for (int i=2; i < 10; i++) {
        massFractions[i-2] = y[i]*mw[i-2];
      }
      massFractions[8] = 0.0;
      GasMixture->setMassFractions(massFractions);
      GasMixture->setState_UV(int_energy, vol, tol);
      y[0] = GasMixture->temperature();
      y[1] = GasMixture->pressure();
    }

    printf("At t = %0.4e      y =%1.16e\n",
           tout, y[0]);

    // Print solutions to file
    for (int i=2; i < 10; i++) {
      std::ofstream outfile;
      outfile.open("Output/test_rk4_" + std::to_string(i-2) + ".txt", std::ios_base::app);
      outfile << std::setprecision(16) << y[i]*mw[i-2] << std::endl;
    }
    std::ofstream outfile;
    outfile.open("Output/test_rk4_temperature.txt", std::ios_base::app);
    outfile << std::setprecision(16) << y[0] << std::endl;
  }

  delete GasMixture;

  return(0);
}
