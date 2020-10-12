/************************************************************************
 * File: test_leap_imex.cpp                                              *
 * Programmers: Cory Mikida                                             *     
 *----------------------------------------------------------------------*
 * Testing Leap IMEX ARK integration of chemistry with PyJac generated source
 * terms, with and without Cantera involvement for dependent variable
 * update.
 * Example problem: homogeneous Cantera mechanism (9-species San-Diego) *
 ************************************************************************/

#include "LeapIMEXMethod.H"
#include "cantera/IdealGasMix.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
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
#define T1    5.0e-9          /* first output time      */
#define TADD  5.0e-9         /* output time factor     */
#define NOUT  19000           /* number of output times */


/***************************** Main Program ******************************/

using namespace std;

main()
{
  double reltol, abstol, tout;
  int iout, flag;
  double y[NSP+1];
  double mw[NSP];
  double massFractions[NSP];
  double corr_norm;
  std::shared_ptr<double> ic(new double [10]);
  std::shared_ptr<double> t_ptr(new double);
  std::shared_ptr<double> dt_ptr(new double);
  dagrt_state_type state;
  double dt_value, t_start;

  // PyJac takes number of moles of each species as input.
  //ic.get()[0] = 1001.0;            /* PyJacV2 convention: temperature goes at the beginning. */
  ic.get()[0] = 1000.0;            /* PyJacV2 convention: temperature goes at the beginning. */
  //ic.get()[1] = 101426.325;            /* PyJacV2 convention: pressure goes at the beginning. */
  ic.get()[1] = 101325.0;            /* PyJacV2 convention: pressure goes at the beginning. */
  ic.get()[2] = 0.0;               /* Leave out last species as well */
  ic.get()[3] = 0.0;
  ic.get()[4] = (0.1949817389200555/0.2072648773462248)/(2*15.9994);
  ic.get()[5] = (0.01228313842616932/0.2072648773462248)/(2*1.00794);
  ic.get()[6] = 0.0;
  ic.get()[7] = 0.0;
  ic.get()[8] = 0.0;
  ic.get()[9] = 0.0;
  
  dt_value = TADD;
  t_start = T0;
  *dt_ptr = dt_value;
  *t_ptr = t_start;

  initialize(state, ic, t_ptr, dt_ptr);

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

  /* Newton iteration using Jacobian, printing results*/
  printf(" \nHomogeneous Cantera problem\n\n");
  Cantera::IdealGasMix * GasMixture;
  //GasMixture = new Cantera::IdealGasMix("Mechanisms/sanDiego.xml");
  GasMixture = new Cantera::IdealGasMix("Mechanisms/SanDiegoN9.xml");
  for (int i=2; i < 10; i++) {
    massFractions[i-2] = ic.get()[i]*mw[i-2];
  }
  massFractions[8] = 0.0;
  GasMixture->setState_TPY(ic.get()[0], ic.get()[1], massFractions);
  double int_energy = GasMixture->intEnergy_mass();
  std::cout << "Internal energy from Cantera of IC: " << std::setprecision(16) << int_energy << std::endl;
  double tol = 1e-10;

  // Remove old copies of output files
  for (int i=2; i < 10; i++) {
    if (remove(("Output/test_leap_" + std::to_string(i-2) + ".txt").c_str()) != 0) {
      perror ("error deleting file");
    }
  }
  if (remove("Output/test_leap_temperature.txt") != 0) {
    perror ("error deleting file");
  }

  /* Timestepping loop */
  for (iout=1, tout=T1; iout <= NOUT; iout++, tout += TADD) {
     
    if (iout == 1) {
      run(state);
      run(state);
    } else {
      run(state);
    } 

    for (int i=2; i < 10; i++) {
      massFractions[i-2] = state.ret_state_y.get()[i]*mw[i-2];
    }
    massFractions[8] = 0.0;
    GasMixture->setMassFractions(massFractions);
    GasMixture->setState_UV(int_energy, vol, tol);
    
    //printf("At t = %0.4e      y =%1.16e\n",
    //       tout, state.ret_state_y.get()[0]);
    printf("At t = %0.4e      y =%1.16e\n",
           tout, GasMixture->temperature());

    // Print solutions to file
    for (int i=2; i < 10; i++) {
      std::ofstream outfile;
      outfile.open("Output/test_leap_" + std::to_string(i-2) + ".txt", std::ios_base::app);
      outfile << std::setprecision(16) << state.ret_state_y.get()[i] << std::endl;
    }
    std::ofstream outfile;
    outfile.open("Output/test_leap_temperature.txt", std::ios_base::app);
    outfile << std::setprecision(16) << state.ret_state_y.get()[0] << std::endl;
  }

  shutdown(state);

  delete GasMixture;

  return(0);
}
