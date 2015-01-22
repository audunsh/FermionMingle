#include <fstream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <string>
#include "integrator/boysfunction.h"
#include "basis/basis.h"
#include "integrator/integrator.h"
#include "solvers/hfsolve.h"
#include "basis/contracted.h"
#include "solvers/rhfsolve.h"
#include "solvers/uhfsolve.h"
#include "solvers/ccsolve.h"
#include "basis/basisbank.h"
#include "interface/fmingle.h"
#include "interface/energyconvert.h"

//double pi = 4*atan(1);


using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
    /**************************************************************************************************************
     * Fermion Mingle Quantum Solver
     * A solver for quantum many-body problems with support for:
     * - Restricted and unrestricted Hartree Fock
     * - Coupled Cluster (CCD, CCSD, CCSDT*)
     * - Gaussian Basis sets
     * - Fermion density evaluation
     *
     * (*) remains to be implemented
     * Written as a university project in C++ by Audun Skau Hansen & Goran Brekke Svaland | Comp-Phys | UiO | 2014
     *
     * Library Armadillo is required for compilation.
     **************************************************************************************************************/
    //Testing UHF procedure
    if(false){
        cout << "Performing an UHF on Li using STO-6G set" << endl;
        vec3 corePos0 = {0,0,0};                            //setting up some position vectors for the cores
        fmingle myparty1;                                   //creating the system
        myparty1.add_nucleus(corePos0, 3);                  //creating a nucleus with charge 1
        myparty1.fminglebasisbank.add_6_311G2df2pd_li(corePos0);  //creating an electrons positioned at core 1 using STO-6G basis set
        myparty1.uhf_solve(1,2);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
    }

    if(true){
        cout << "Testing the electrongas" << endl;
        fmingle myparty;
        myparty.add_electrongas(1,1);
        myparty.BS.init_electron_integrals();
        energyconvert energy2(myparty.BS.gasbasis.eref(14));
        energyconvert energy_a(myparty.BS.gasbasis.analytic_energy(1));
        //cout << "Reference energy:" << energy2.as("rydberg")/14.0<< endl;

        //myparty.BS.h.print();
        //myparty.BS.v.print();
        //myparty.printing = true;
        //myparty.rhf_solve(4);
        //myparty.fminglesolver_hf.solve_rhf(2.0);
        //myparty.fminglesolver_rhf.coupledMatrix.print();
        //cout << myparty.fminglesolver_rhf.energyCalc()<< endl;

        //rhfsolve ecalc (myparty.BS, 14);


        cout << "Calculating energy" << endl;
        //energyconvert energy(ecalc.energyCalc());

        cout << "Analytic energy    :" << energy_a.rydberg() << endl;
        cout << "Resulting HF energy:" << energy2.rydberg()/14.0 << endl; //in Rydbergs/electron
        cout << "Compared to        :" << "1.93434 rydbergs per particle from G.Baardsen." << endl;
        cout << "With a ratio of    :" << (energy2.rydberg()/14.0)/1.93434 << endl;

    }

    if(false){
        cout << "Performing UHF+CCSD on H2 using STO-6G set" << endl;
        vec3 corePos1 = {0,0,0};                            //setting up some position vectors for the cores
        vec3 corePos2 = {0,0,6.0};
        fmingle myparty2;                                   //creating the system
        myparty2.printing = true;
        myparty2.add_nucleus(corePos1, 1);                  //creating a nucleus with charge 1
        myparty2.add_nucleus(corePos2, 1);                //creating a nucleus with charge 1
        myparty2.fminglebasisbank.add_STO_3G_h(corePos1);  //creating an electrons positioned at core 1 using STO-6G basis set
        myparty2.fminglebasisbank.add_STO_3G_h(corePos2); //creating an electrons positioned at core 2 using STO-6G basis set
        //myparty2.rhf_solve(2);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons

        myparty2.uhf_solve(1,1);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
        //myparty2.ccsd_solve(2);

    }

    if(false){
        cout << "Performing RHF+CCSD on H2 using STO-6G set" << endl;
        vec3 corePos1 = {0,0,0};                            //setting up some position vectors for the cores
        vec3 corePos2 = {0,0,2.287};

        fmingle myparty2;                                   //creating the system
        myparty2.printing = true;
        myparty2.add_nucleus(corePos1,8);                  //creating a nucleus with charge 1
        //myparty2.add_nucleus(corePos2,8);                //creating a nucleus with charge 1
        myparty2.fminglebasisbank.add_STO_3G_o(corePos1);  //creating an electrons positioned at core 1 using STO-6G basis set
        //myparty2.fminglebasisbank.add_STO_3G_o(corePos2); //creating an electrons positioned at core 2 using STO-6G basis set
        myparty2.rhf_solve(8);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
        myparty2.ccsd_solve(8);

        //cout << "Number of primitives:" << myparty2.fminglebasisbank.bs.nPrimitivesInState(0) << endl;

        /*
        myparty2.reset();
        cout << endl << "Performing UHF on H2 using STO-6G set" << endl;
        //cout << "Number of primitives:" << myparty2.fminglebasisbank.bs.nPrimitivesInState(0) << endl;
        myparty2 = fmingle();
        myparty2.printing = true;
        myparty2.add_nucleus(corePos1, 1);                  //creating a nucleus with charge 1
        myparty2.add_nucleus(corePos2, 1);                //creating a nucleus with charge 1
        myparty2.fminglebasisbank.add_STO_6G_h(corePos1);  //creating an electron positioned at core 1 using STO-6G basis set
        myparty2.fminglebasisbank.add_STO_6G_h(corePos2); //creating an electron positioned at core 2 using STO-6G basis set
        myparty2.uhf_solve(1,1);                            //perform a unrestricted hartree fock procedure for 1 up electron, 2 down electrons
        //myparty2.ccsd_solve(2);

        */
        cout << "Correlation energy:" << myparty2.correlation_energy << endl;
        cout << "RHF Energy:" << myparty2.rhf_energy << endl;

    }

    if(false){
        double x0 = 0.1;
        double x1 = 6.0;
        int N = 50;
        double dx = (x1-x0)/(N-1);
        mat energyVsPosition = zeros(N,2);
        //vec enrg = zeros(N);
        vec3 corePos1 = {0,0,0};
        vec3 corePos2 = {0,0,0};

        fmingle myparty2;                                   //creating the system

        for(int i = 0; i<N ; i++){
            myparty2 = fmingle ();

            energyVsPosition(i,0) = x0 + dx*i;
            corePos1 = {0,0,0};
            corePos2 = {0,0,energyVsPosition(i,0)};

            myparty2.fminglebasisbank.bs.add_nucleus(corePos1, 1);
            myparty2.fminglebasisbank.bs.add_nucleus(corePos2, 1);
            myparty2.fminglebasisbank.add_STO_6G_h(corePos1);  //creating an electron positioned at core 1 using STO-6G basis set
            myparty2.fminglebasisbank.add_STO_6G_h(corePos2);  //creating an electron positioned at core 2 using STO-6G basis set

            myparty2.rhf_solve(2);                            //perform a restricted hartree fock procedure for 1 up electron, 2 down electrons
            myparty2.ccsd_solve(2);

            //energyVsPosition(i,1) = myparty2.rhf_energy + myparty2.correlation_energy;
            energyVsPosition(i,1) = myparty2.rhf_energy + myparty2.correlation_energy;

            //energyVsPosition(i,1) = myparty2.fminglesolver_cc.t2c(3,2)(1,0);

            //cout << myparty2.fminglesolver_hf. << endl;
            if(i>1){
                if(abs(energyVsPosition(i-1,1) - energyVsPosition(i,1))>1.0){
                    cout << "A convergence likely occured." << endl;
                    //cout << myparty2.fminglesolver_rhf.iterations << endl;
                }
            }
            myparty2.reset();
            //myparty.rhf_solve(2);
            //e(i) = myparty.rhf_energy;

        }
        energyVsPosition.print();
        energyVsPosition.save("STO6G_h2_RHF_CCSD_N40_t1_.txt", raw_ascii);
    }

    return 0;

}
