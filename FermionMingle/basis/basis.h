#ifndef BASIS_H
#define BASIS_H

#include <string>
#include <armadillo>
#include "basis/contracted.h"
#include "basis/primitive.h"
#include "basis/electrongas.h"
#include "integrator/integrator.h"
#include "integrator/boysfunction.h"


using namespace std;
using namespace arma;

class basis
{
public:
    basis();
    void set_size(int N);                                         //set size of basis, N=number of orbitals
    void reset(); //reset all basis information
    void read(string filename, int Zn);                           //read basis from file
    void set_orthonormal();                                       //if true: set the overlap matrix to I
    void expand();                                                //expand basis for explicit spin-dependence
    void init_integrals();                                        //set up and solve all integrals for a gaussian basis
    void init_STO_3G(string configuration, double nProtons);                       //initialize STO-3G basis set for given configuration
    void init_HTO4(double nProtons);                                //initialize hydrogen type basis with 4 orbitals
    void init_molecule(string configuration, vec nProtons, field<vec> corePos);
    void printAllContracted();
    void flip_index(); //Change indexing to match Thijssen

    //Electron Gas Implementation
    void setup_electrongas(int NEv, double Lv);
    void init_electron_integrals();
    int L, NE;
    electrongas gasbasis;

    //implementing a more streamlined interface for dealing with the basis
    void add_state();
    void add_primitive_to_state(int Stateindex, Primitive P);
    void add_nucleus(vec3 pos, int charge);
    void add_STO_NG_s_orbital(int NG, vec weigths, vec exponents, vec3 corePos);
    void add_STO_NG_p_orbital(int NG, vec exponents, vec weigths, vec3 corePos);
    void add_STO_NG_d_orbital(int NG, vec exponents, vec weigths, vec3 corePos);
    //void import(TurboMoleParser system, vec3 corePos);

    void init_H2(vec3 corePos1, vec3 corePos2);
    void init_Be2(vec3 corePos1, vec3 corePos2);
    void init_H2O(vec3 corePos1, vec3 corePos2, vec3 corePos3);
    void init_H2Be(vec3 corePos1, vec3 corePos2, vec3 corePos3);
    void init_O2(vec3 corePos1, vec3 corePos2);
    void init_O(vec3 corePos1);
    void init_Ne(vec3 corePos1);
    void init_He(vec3 corePos1);
    void init_Be(vec3 corePos1);
    void init_He2(vec3 corePos1, vec3 corePos2);



    void add_atom_STO3G(string configuration, vec3 corePos);
    //void init_Be2(vec3 corePos1, vec3 corePos2);

    double state(int p, int q, int r, int s, double D, double E); //function to evaluate spin-dependence
    int Nstates, Nstates2;                                        //number of basis states (spin not included)
    field<mat> v;                                                 //basis, spin not included
    field<mat> V;                                                 //basis with spin included
    double Z;                                                     //interaction parameter (taken to be number of protons)
    double h0(int i, int j);                                      //one body interaction
    mat S,h,H,nuclearPotential;                                   //overlap matrix, onebody interaction matrices
    double nnInteraction(); //nucleon-nucleon contribution to energy

    //Correction due to normalization factor from the Turbomole format
    Primitive turbomolePrimitive(double weight, double exponent, double i, double j, double k, vec3 corePos);
    double factorial(double n);
    double turboNormalization(double x, double i, double j ,double k);
    double evaluateContracted(int n, vec3 r);

    //Need to access this externally, temporarily
    vec nucleusCharges;
    vec nPrimitivesInState;


private:
    contracted basisSet[]; //use one of these...
    vector<contracted> basisSts; //use one of these...
    field<vec3> nucleusPositions;

    int nNucelons;
    int Nprimitives;
    double pi = 4*atan(1);

};

#endif // BASIS_H
