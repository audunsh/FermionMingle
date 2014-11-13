#ifndef HFSOLVE_H
#define HFSOLVE_H
#include <armadillo>
#include <string>
#include "basis/basis.h"
#include "solvers/uhfsolve.h"
#include "solvers/rhfsolve.h"

using namespace std;
using namespace arma;

class HFSolve{
    int Z,N, Nstates;
public:
    HFSolve();
    HFSolve(basis BS);
    double solve(int N_electrons);
    void solve_rhf(int N_electrons);
    void solve_uhf(int N_electrons_up, int N_electrons_down);
    void reset();
    //rhfsolve restrictedhf;
    //uhfsolve urestrictedhf;
    basis Bs;


    //The following variables must be available for the CCSolve class when performing coupled cluster calculations
    mat C;          //Coefficient matrix
    mat Cu;          //Coefficient matrix
    mat Cd;          //Coefficient matrix
    mat F;          //Fock matrix
    vec epsilon;    //eigenvalues from current diagonalization
    mat P;          //Density matrix
    mat Cprime; //transformed Coefficient matric
    double energy;

private:

};

#endif // HFSOLVE_H
