#include "interface/fmingle.h"
#include <iostream>
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

fmingle::fmingle()
{
    //this class will be the main user interface for the solver
    //A typical runtime should only communicate with the respective solvers through this class
    basis BS;
    initialized = 0;
    printing = false;
    basisbank wrapped (BS);
    fminglebasisbank = wrapped;
    fminglebasisbank.bs.Nstates = 0;
    report = "";
}

void fmingle::add_nucleus(vec3 corePos, int charge){
    fminglebasisbank.bs.add_nucleus(corePos, charge);
}

void fmingle::add_orbitals(vec3 corePos, string config){
}

void fmingle::initialize(){
    fminglebasisbank.bs.set_size(fminglebasisbank.bs.Nstates);
    fminglebasisbank.bs.init_integrals();
    HFSolve solverobject(fminglebasisbank.bs);
    fminglesolver_hf = solverobject;
}

void fmingle::rhf_solve(int nElectrons){
    if(initialized == 0){initialize();
    }
    initialized =1;
    fminglesolver_hf.solve_rhf(nElectrons);
    rhf_energy = fminglesolver_hf.energy;
    if(printing){
    cout << "-------------------------------------------------------------------" << endl;
    cout << std::setprecision(14)<<"  Restricted Hartree-Fock energy:" << rhf_energy << endl;
    cout << "-------------------------------------------------------------------" << endl;}
};

void fmingle::uhf_solve(int nElectronsUp, int nElectronsDown){
    if(initialized == 0){initialize();
    }
    initialized =2;
    fminglesolver_hf.solve_uhf(nElectronsUp, nElectronsDown);
    uhf_energy = fminglesolver_hf.energy;
    if(printing){
    cout << "-------------------------------------------------------------------" << endl;
    cout << std::setprecision(14)<<"Unrestricted Hartree-Fock energy:" << uhf_energy << endl;
    cout << "-------------------------------------------------------------------" << endl;}
};

void fmingle::ccsd_solve(int nElectrons){
    if(initialized==0){
        cout << "The basis is not initialized." << endl;
        //Do nothing.
    }
    if(initialized==1){
        //Perform ccd for a RHF basis
        ccsolve solverobject(fminglesolver_hf);
        fminglesolver_cc = solverobject;
        fminglesolver_cc.init_RHF_basis();
        fminglesolver_cc.CCSD(nElectrons);
        if(printing){
        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<< "       CCSD Electron correlation:" << fminglesolver_cc.correlation_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;

        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<<"                    Total energy:" << fminglesolver_cc.correlation_energy+rhf_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;}
        correlation_energy = fminglesolver_cc.correlation_energy;


    }
    if(initialized==2){
        //Perform ccd for a RHF basis
        //cout << "Entering UHF CCSD procedure."<< endl;
        ccsolve solverobject(fminglesolver_hf);
        fminglesolver_cc = solverobject;
        fminglesolver_cc.init_UHF_basis();
        fminglesolver_cc.CCSD(nElectrons);
        correlation_energy = fminglesolver_cc.correlation_energy;
        if(printing){
        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<<"       CCSD Electron correlation:" << fminglesolver_cc.correlation_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;
        cout << "-------------------------------------------------------------------" << endl;
        cout << std::setprecision(14)<<"                    Total energy:" << fminglesolver_cc.correlation_energy+uhf_energy << endl;
        cout << "-------------------------------------------------------------------" << endl;}
    }


}

void fmingle::reset(){
    //fminglebasisbank.bs = basis ();
    //fminglebasisbank.bs.set_size(0);
    fminglebasisbank.bs.reset();
    fminglebasisbank.bs.set_size(0);
    initialized = 0;
}

void fmingle::sweep_h2(double x0, double x1, int N){
    //perform a sweep for h2
    double dx = (x1-x0)/(N-1);
    vec xpos = zeros(N);
    vec enrg = zeros(N);
    vec3 corePos1 = {0,0,0};
    vec3 corePos2 = {0,0,0};
    for(int i = 0; i<N ; i++){
        //wrapped = basisbank();
        //BS = basis();
        //basisbank wrapped (BS);

        cout << "Reset basisbank" << endl;
        //fminglebasisbank.bs = basis();
        fminglebasisbank.bs.reset();
        xpos(i) = x0 + dx*i;
        corePos1 = {0,0,0};
        corePos2 = {0,0,xpos(i)};

        fminglebasisbank.bs.add_nucleus(corePos1, 1);
        fminglebasisbank.bs.add_nucleus(corePos2, 1);
        fminglebasisbank.add_STO_6G_h(corePos1);  //creating an electron positioned at core 1 using STO-6G basis set
        fminglebasisbank.add_STO_6G_h(corePos2);  //creating an electron positioned at core 2 using STO-6G basis set
        initialize();
        //myparty.rhf_solve(2);
        //e(i) = myparty.rhf_energy;

    }


}

void fmingle::ccd_solve(int nElectrons){

}
