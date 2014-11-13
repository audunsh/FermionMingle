#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "basis/primitive.h"
#include "integrator/boysfunction.h"

using namespace std;
using namespace arma;

class integrator
{
public:
    integrator(Primitive &pA, Primitive &pB, BoysFunction &boysf);
    double overlapIntegral(Primitive &pA, Primitive &pB);

    //the framework for the recursive algorithms
    //void setupHermiteCoefficients();
    void setupEij();
    void setupEcd();
    void setupRtuv(vec3 &nucleiPos);
    void setupRtau();


    //the integrals
    double overlap();
    double kinetic();
    double pNuclei();
    double pp(Primitive &pC, Primitive &pD);

private:
    //double pi = 4*atan(1);
    field <cube> Eij;
    field <cube> Ecd;
    field <cube> Rtuv;
    field <cube> Rtau;



    vec3 P, pAijk, pBijk, pCijk, pDijk,A,B,C,D,Xab,Xcd,Xpa,Xpb,Xqc,Xqd,Rpc,Sijk,Tijk,Q, Rpq;
    double a,b,c,d,p,mu, Xab2,Xcd2,wA,wB,wC,wD,R, Rpc2,Rpq2,S,q, alpha;
    BoysFunction boys;
};

#endif // INTEGRATOR_H
