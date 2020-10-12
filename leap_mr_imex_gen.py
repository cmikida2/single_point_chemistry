from __future__ import print_function

import dagrt.codegen.cxx as cxx
from leap.multistep.multirate import (MultiRateMultiStepMethodBuilder,
                                      MultiRateHistory as MRHistory,
                                      rhs_policy)
from time_int_helpers import am_solver_hook


def main():
    # order = 4
    # hist_length = 4
    order = 3
    hist_length = 3
    static_dt = True
    stepper = MultiRateMultiStepMethodBuilder(
                order,
                (
                    (
                        'dt', 'fast', '=',
                        MRHistory(1, "<func>f", ("fast", "slow"),
                                  hist_length=hist_length,
                                  is_rhs_implicit=True),
                        ),
                    (
                        'dt', 'slow', '=',
                        MRHistory(1, "<func>s", ("fast", "slow"),
                                  rhs_policy=rhs_policy.late,
                                  hist_length=hist_length,
                                  is_rhs_implicit=False),
                        ),),
                static_dt=static_dt)

    from dagrt.function_registry import (
            base_function_registry, register_function, UserType)

    # This stays unchanged from the existing explicit RK4 RHS.
    freg = register_function(base_function_registry, "<func>s",
                             ("t", "fast", "slow"), result_names=("result",),
                             result_kinds=(UserType("slow"),))
    freg = freg.register_codegen("<func>s", "cxx",
                                 cxx.CallCode("""

                // Purely homogeneous chemistry simulation.
                for (int i = 0;i < 4;i++){

                 ${result}[i] = 0.0;

                }

                """))

    # Here, we need to call a new chemistry RHS that loops through
    # all the points *under the hood.*
    freg = register_function(freg, "<func>f", ("t", "fast", "slow"),
                             result_names=("result",),
                             result_kinds=(UserType("fast"),))
    freg = freg.register_codegen("<func>f", "cxx",
                                 cxx.CallCode("""

                // PyJac inputs.
                double jac[(NS+1)*(NS+1)];
                double jac_trans[(NS+1)*(NS+1)];
                double phi[(NS+1)];
                double phi_guess[(NS+1)];
                double phi_old[(NS+1)];
                double dphi[(NS+1)];
                double dphi_old[(NS+1)];
                double corr[(NS+1)];
                double corr_weights[(NS+1)];
                double corr_weighted[(NS+1)];
                double reltol = 1e-6;
                double abstol = 1e-12;
                /* For Lapack */
                int ipiv[NS+1], info;
                int nrhs = 1;
                int nsp_l = NS+1;
                // Dummy time for pyJac.
                double tout = 0;

                // Work array for PyJac.
                double* rwk_dphi = (double*)malloc(245 * sizeof(double));
                memset(rwk_dphi, 0, 245 * sizeof(double));
                double massFractions[NS];
                double mw[NS];

                // FIXME: Assumes the last (inert) species is nitrogen.
                mw[NS-1] = 2*14.00674;
                for (int i = 0; i < NS-1; ++i ){
                  mw[i] = mw[NS-1]*mw_factor[i];
                }

                double rho = 0.2072648773462248;
                double vol = 1.0 / rho;
                double tol = 1e-10;
                double mass_sum = 0.0;
                for (int i = 2; i <= NS; ++i ){
                  massFractions[i-2] = ${fast}[i]*mw[i-2];
                  mass_sum += massFractions[i-2];
                }
                //massFractions[8] = 1 - mass_sum;
                massFractions[8] = 0.0;

                // PyJac converted input state.
                //phi[0] = ${fast}[0];
                //phi[1] = ${fast}[1];
                // Update temperature and pressure using Cantera?
                Cantera::IdealGasMix * GasMixture;
                GasMixture = new Cantera::IdealGasMix("Mechanisms/sanDiego.xml");
                double int_energy = 788261.179011143;
                GasMixture->setMassFractions(massFractions);
                GasMixture->setState_UV(int_energy, vol, tol);
                phi[0] = GasMixture->temperature();
                phi[1] = GasMixture->pressure();
                delete GasMixture;
                for (int j=2;j <= NS; j++) {
                  phi[j] = massFractions[j-2]/mw[j-2];
                }
                // PyJac call for source term
                species_rates (&tout, &vol, phi, dphi, rwk_dphi);

                for (int j=0;j <= NS; j++) {
                  ${result}[j] = dphi[j];
                }

                """))

    # The tricky part - this is going to require a
    # nonlinear solve of some kind.
    freg = register_function(freg, "<func>solver",
                             ("fast", "slow", "coeff", "t"),
                             result_names=("result",),
                             result_kinds=(UserType("fast"),))
    freg = freg.register_codegen("<func>solver", "cxx",
                                 cxx.CallCode("""

                // PyJac inputs.
                double jac[(NS+1)*(NS+1)];
                double jac_trans[(NS+1)*(NS+1)];
                double jac_sub[(NS-1)*(NS-1)];
                double phi[(NS+1)];
                double phi_guess[(NS+1)];
                double phi_old[(NS+1)];
                double dphi[(NS+1)];
                double dphi_sub[(NS-1)];
                double dphi_old[(NS+1)];
                double corr[(NS+1)];
                double corr_weights[(NS+1)];
                double corr_weighted[(NS+1)];
                double reltol = 1e-6;
                double abstol = 1e-12;
                /* For Lapack */
                int ipiv[NS-1], info;
                int nrhs = 1;
                int nsp_l = NS-1;
                int nsp_l_full = NS+1;
                // Dummy time for pyJac.
                double tout = 0;

                // Work array for PyJac.
                double* rwk_dphi = (double*)malloc(245 * sizeof(double));
                memset(rwk_dphi, 0, 245 * sizeof(double));
                double* rwk_jac = (double*)malloc(245 * sizeof(double));
                memset(rwk_jac, 0, 245 * sizeof(double));

                /* 1D point-sized identity matrix */
                double ident[(NS-1)*(NS-1)];
                for (int i=0; i < (NS-1)*(NS-1); i++) {
                    if (i % (NS) == 0) {
                        ident[i] = 1;
                    } else {
                        ident[i] = 0;
                    }
                }

                // THIS IS WHERE ALL OF THE
                // NEWTON/PYJAC STUFF GOES.
                // Get chemical state at this point.

                double massFractions[NS];
                double mw[NS];

                // FIXME: Assumes the last (inert) species is nitrogen.
                mw[NS-1] = 2*14.00674;
                for (int i = 0; i < NS-1; ++i ){
                  mw[i] = mw[NS-1]*mw_factor[i];
                }

                double rho = 0.2072648773462248;
                double mass_sum = 0.0;
                for (int i = 2; i <= NS; ++i ){
                  massFractions[i-2] = ${fast}[i]*mw[i-2];
                  mass_sum += massFractions[i-2];
                }
                //massFractions[8] = 1 - mass_sum;
                massFractions[8] = 0.0;

                double vol = 1.0 / 0.2072648773462248;
                double tol = 1e-10;

                // PyJac converted input state.
                //phi[0] = ${fast}[0];
                //phi[1] = ${fast}[1];
                Cantera::IdealGasMix * GasMixture;
                GasMixture = new Cantera::IdealGasMix("Mechanisms/sanDiego.xml");
                double int_energy = 788261.179011143;
                GasMixture->setMassFractions(massFractions);
                GasMixture->setState_UV(int_energy, vol, tol);
                phi[0] = GasMixture->temperature();
                phi[1] = GasMixture->pressure();
                //delete GasMixture;
                for (int j=2;j <= NS; j++) {
                  phi[j] = massFractions[j-2]/mw[j-2];
                }
                for (int j=0;j <= NS; j++) {
                  phi_old[j] = phi[j];
                }
                for (int j=0;j <= NS; j++) {
                  phi_guess[j] = phi[j];
                }
                // Newton loop within this point (for now).
                double corr_norm = 1.0;
                // PyJac call for source term
                species_rates (&tout, &vol, phi, dphi_old, rwk_dphi);
                //while (abs(corr_norm) >= reltol) {
                while (abs(corr_norm) >= 1e-8) {
                  species_rates (&tout, &vol, phi_guess, dphi, rwk_dphi);
                  // Get Jacobian at this point.
                  jacobian (&tout, &vol, phi_guess, jac, rwk_jac);
                  for (int j=0;j <= NS; j++) {
                    // IMEX AM
                    dphi[j] = phi_guess[j] - phi_old[j] - ${coeff} * dphi[j];
                  }

                  // Take subset of RHS vector.
                  for (int j=2;j <= NS; j++) {
                    dphi_sub[j-2] = dphi[j];
                  }

                  // Transpose the Jacobian.
                  // PYJAC outputs Fortran-ordering
                  for (int i = 0; i < (NS+1); ++i )
                  {
                    for (int j = 0; j < (NS+1); ++j )
                    {
                      // Index in the original matrix.
                      int index1 = i*(NS+1)+j;

                      // Index in the transpose matrix.
                      int index2 = j*(NS+1)+i;

                      jac_trans[index2] = jac[index1];
                    }
                  }

                  for (int i=0; i<(NS+1)*(NS+1); i++) {
                    jac[i] = jac_trans[i];
                  }

                  // Take subset of Jacobian for algebraic changes.
                  for (int i = 0; i < (NS-1); ++i )
                  {
                    for (int j = 0; j < (NS-1); ++j )
                    {
                      jac_sub[i*(NS-1)+j] = jac[(i+2)*(NS+1)+2+j];
                    }
                  }

                  // Make the algebraic changes,
                  // using PyJac routines as needed.
                  // Get internal energies.
                  double int_energies[NS];
                  eval_u(phi_guess, int_energies);
                  // Get CVs.
                  double cvs[NS];
                  eval_cv(phi_guess, cvs);
                  // Get total CV.
                  double mass_sum = 0.0;
                  for (int i = 2; i <= NS; ++i ){
                    massFractions[i-2] = phi_guess[i]*mw[i-2];
                    mass_sum += massFractions[i-2];
                  }
                  //massFractions[8] = 1 - mass_sum;
                  massFractions[8] = 0.0;
                  double cv_total = 0.0;
                  for (int i = 0; i <= NS-1; ++i ){
                    cv_total += cvs[i]*(massFractions[i]/mw[i]);
                  }

                  // Modify the Jacobian for the algebraic constraint.
                  for (int i = 0; i < (NS-1); ++i )
                  {
                    for (int j = 0; j < (NS-1); ++j )
                    {
                      jac_sub[i*(NS-1)+j] -= jac[j+2]*int_energies[i]/cv_total;
                    }
                  }
                  /* Subtract from identity to get jac for Newton */
                  for (int j=0;j < (NS-1)*(NS-1); j++) {
                    jac_sub[j] = ident[j] - ${coeff} * jac_sub[j]; // IMEX AM
                    jac_sub[j] = -jac_sub[j];
                  }
                  // Do the inversion with the assistance of Lapack.
                  dgesv_(&nsp_l, &nrhs, jac_sub, &nsp_l, ipiv, dphi_sub, &nsp_l, &info);
                  // Add the correction to the buffer.
                  mass_sum = 0.0;
                  for (int i = 2; i <= NS; ++i ){
                    massFractions[i-2] = (phi_guess[i] + dphi_sub[i-2])*mw[i-2];
                    mass_sum += massFractions[i-2];
                  }
                  //massFractions[8] = 1 - mass_sum;
                  massFractions[8] = 0.0;
                  GasMixture->setMassFractions(massFractions);
                  GasMixture->setState_UV(int_energy, vol, tol);
                  corr[0] = GasMixture->temperature() - phi_guess[0];
                  corr[1] = GasMixture->pressure() - phi_guess[1];
                  for (int j=2;j <= NS; j++) {
                    corr[j] = dphi_sub[j-2];
                  }
                  for (int j=0;j <= NS; j++) {
                    corr_weights[j] = 1.0 / (reltol * abs(phi_guess[j]) + abstol);
                    corr_weighted[j] = corr[j]*corr_weights[j];
                  }

                  //corr_norm = dnrm2_(&nsp_l, corr, &nrhs);
                  corr_norm = dnrm2_(&nsp_l_full, corr_weighted, &nrhs);
                  //std::cout << "Correction norm: " << corr_norm << std::endl;
                  for (int j=0;j <= NS; j++) {
                    phi_guess[j] = phi_guess[j] + corr[j];
                  }
                }
                // Now outside the Newton loop (presumably having converged),
                // we update the RHS using the actual state.
                species_rates (&tout, &vol, phi_guess, dphi, rwk_dphi);
                for (int j=0;j <= NS; j++) {
                  ${result}[j] = dphi[j];
                }
                delete GasMixture;

                """))

    code = stepper.generate()

    # Implicit solve thingy
    from leap.implicit import replace_AssignImplicit
    code = replace_AssignImplicit(code, {"solve": am_solver_hook})

    # print(code)
    codegen = cxx.CodeGenerator(
            'LeapIMEX',
            user_type_map={
                "fast": cxx.ArrayType((10,), cxx.BuiltinType('double'),),
                "slow": cxx.ArrayType((4,), cxx.BuiltinType('double'),),
                },
            function_registry=freg,
            emit_instrumentation=True,
            timing_function="clock",
            header_preamble="\n#include \"mechanism.hpp\"\n#include " +
                            "\"species_rates.hpp\"\n#include " +
                            "\"jacobian.hpp\"\n#include \"memcpy_2d.hpp" +
                            "\"\n#include \"lapack_kernels.H\"\n#include " +
                            "\"cantera/IdealGasMix.h\"\n#include " +
                            "\"cantera/thermo.h\"\n#include " +
                            "\"cantera/kinetics.h\"\n#include " +
                            "\"cantera/transport.h\"")

    import sys

    # Write out Leap/Dagrt code:
    with open(sys.argv[1], "a") as outf:
        code_str = codegen(code)
        print(code_str, file=outf)


if __name__ == "__main__":
    main()
