#ifndef UHFSOLVE_H
#define UHFSOLVE_H

#include "basis/basis.h"

using namespace std;
using namespace arma;

class uhfsolve
{
public:
    uhfsolve();
    uhfsolve(basis BS, int N_spin_up, int N_spin_down);
    void setupUnitMatrices(); //Bring overlap matrix to unit form
    void setupP();       //setup density matrix, make a first guess
    void setupF();       //setup the Fock matrix
    void diagonalizeF(); //diagonalize the Fock matrix
    void normalizeC();   //normalize the coefficient matrix
    void updateP();      //construct new density matrix
    void setupCoupledMatrix(); //set up direct and exchange terms from basis
    void printMatrices();//print the iterating matrices
    void setupTotalCoefficientMatrix();
    bool convergenceCriteria(); //check for convergence
    double energy();     //return ground state energy
    double solve();      //automated solving process, returns energy
    double coupledMatrixTilde(int p, int q, int r, int s); //return direct and exchange term

    //new code to debug, possibly for deletion
    void setupCoupledMatrix_unused();
    double calcEnergy2();
    double energyCalc();
    double evaluateProbabilityDensity(vec3 r);
    void createDensityMap(string filename);
    void reset(basis BS, int N, int Z);

    //Trying to make some objects availabe to external classes, rammifications unknown
    mat C;          //Coefficient matric
    mat Cu;          //Coefficient matrix u
    mat Cd;          //Coefficient matrix d
    field<mat> coupledMatrix;
    mat F;          //Fock matrix
    mat Fu;          //Fock matrix up
    mat Fd;          //Fock matrix down
    int nElectrons; //number of electrons in total
    int nElectronsU; //number of electrons with spin-up config
    int nElectronsD; // --" -- spin-down
    int nStates;    //number of states
    basis Bs;
    vec epsilon;    //eigenvalues from current diagonalization

    vec epsilonu;    //eigenvalues from current diagonalization
    vec epsilond;    //eigenvalues from current diagonalization

    mat P;          //Density matrix
    mat Pu;          //Density matrix
    mat Pd;          //Density matrix

    mat Cprime; //transformed Coefficient matric
    mat Cprimeu; //transformed Coefficient matric
    mat Cprimed; //transformed Coefficient matric


private:
    cube densityMap;

    mat V;      //Transformation matrix

    mat P_prev; //previous density matrix
    mat Pu_prev; //previous density matrix
    mat Pd_prev; //previous density matrix


    mat Fprime; //transformed Fock matrix
    mat Fprimeu; //transformed Fock matrix
    mat Fprimed; //transformed Fock matrix


    mat G;      //Coulomb and exchange contribution matrix
    mat U;      //unitary matrix


    vec epsilon_prev; //eigenvalues from previous diagonalization
    vec epsilonu_prev; //eigenvalues from previous diagonalization
    vec epsilond_prev; //eigenvalues from previous diagonalization
    vec s_diag;       //diagonalized S matrix
    int iterations; //number of iterations (counter)


    int nProtons;   //number of protons, can be removed


    double energyPrev = 10e10;
    double tolerance = 10e-10; //-10
    double dampingFactor = 0.99; //0.95
};

#endif // UHFSOLVE_H
