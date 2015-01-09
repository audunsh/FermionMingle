#ifndef RHFSOLVE_H
#define RHFSOLVE_H

#include "basis/basis.h"
#include <armadillo>

class rhfsolve
{
public:
    rhfsolve();
    rhfsolve(basis BS, int N);
    void setupUnitMatrices(); //Bring overlap matrix to unit form
    void setupP();       //setup density matrix, make a first guess
    void setupF();       //setup the Fock matrix
    void diagonalizeF(); //diagonalize the Fock matrix
    void normalizeC();   //normalize the coefficient matrix
    void updateP();      //construct new density matrix
    void setupCoupledMatrix(); //set up direct and exchange terms from basis
    void printMatrices();//print the iterating matrices
    bool convergenceCriteria(); //check for convergence
    double energy();     //return ground state energy
    double solve();      //automated solving process, returns energy
    double coupledMatrixTilde(int p, int q, int r, int s); //return direct and exchange term
    double getOrbitalEnergy(int i);

    //new code to debug, possibly for deletion
    void setupCoupledMatrix_unused();
    double calcEnergy2();
    double energyCalc();
    double evaluateProbabilityDensity(vec3 r);
    void createDensityMap(string filename);
    void reset(basis BS, int N);

    //Trying to make some objects availabe to external classes, rammifications unknown
    mat C;          //Coefficient matric
    field<mat> coupledMatrix;
    mat F;          //Fock matrix
    int nElectrons; //number of electrons
    int nStates;    //number of states
    basis Bs;
    vec epsilon;    //eigenvalues from current diagonalization
    mat P;          //Density matrix
    mat Cprime; //transformed Coefficient matric
    int iterations; //number of iterations (counter)

private:
    cube densityMap;

    mat V;      //Transformation matrix

    mat P_prev; //previous density matrix


    mat Fprime; //transformed Fock matrix


    mat G;      //Coulomb and exchange contribution matrix
    mat U;      //unitary matrix


    vec epsilon_prev; //eigenvalues from previous diagonalization
    vec s_diag;       //diagonalized S matrix


    int nProtons;   //number of protons, can be removed


    double energyPrev = 10e10;
    double tolerance = 10e-10; //-10
    double dampingFactor = 0.95; //0.95
};

#endif // rhfsolve_H
