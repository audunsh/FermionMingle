#include "integrator/integrator.h"
#include "basis/primitive.h"
#include "integrator/boysfunction.h"
#include <armadillo>

using namespace std;
using namespace arma;

double pi = 4*atan(1);

/***********************************************************************************
 *  Calculate integrals between primitive objects collected from class Primitive.
 *
 **********************************************************************************/
integrator::integrator(Primitive &pA, Primitive &pB, BoysFunction &boysf){
    boys = boysf;

    a = pA.exponent();          // exponential constant.
    A = pA.nucleusPosition();   // nucleus pA position
    pAijk(0) = pA.xExponent();
    pAijk(1) = pA.yExponent();
    pAijk(2) = pA.zExponent();
    wA = pA.weight();

    b = pB.exponent();          // exponential constant.
    B = pB.nucleusPosition();   // nucleus pB position
    pBijk(0) = pB.xExponent();
    pBijk(1) = pB.yExponent();
    pBijk(2) = pB.zExponent();
    wB = pB.weight();

    //Setting up shared variables
    p = a+b;
    P = (a*A+b*B)/p;
    Xab = A-B;
    Xpa = P-A;
    Xpb = P-B;
    Xab2 = Xab(0)*Xab(0)+Xab(1)*Xab(1)+Xab(2)*Xab(2);

    setupEij();
}


void integrator::setupEij(){
    Eij.set_size(3);
    int I,J,T;

    for(int coord=0;coord<3;coord++){
        I = pAijk(coord);
        J = pBijk(coord);
        T = pAijk(coord)+pBijk(coord);
        Eij(coord).set_size(I+5,J+5,T+6); //not sure about these values (could be lower)
        Eij(coord).zeros();
        Eij(coord) (1,1,1) = exp(-(a*b/p)*(Xab(coord)*Xab(coord)));
        //Eij(coord) (1,1,1) = exp(-(a*b/p)*Xab2);
        //cout << exp(-(a*b/p)*Xab2) << endl;
        for(int i=1;i<I+4;i++){
            for(int t=1;t<T+5; t++){
                Eij(coord) (i+1,1,t) = Eij(coord) (i,1,t-1)/(2*p) + Xpa(coord)*Eij(coord) (i,1,t) + t* Eij (coord) (i,1,t+1);
            }
        }
        for(int j=1;j<J+4;j++){
            for(int i=1;i<I+4;i++){
                for(int t=1;t<T+5; t++){
                    Eij(coord) (i,j+1,t) = Eij(coord) (i,j,t-1)/(2*p) + Xpb(coord)*Eij(coord) (i,j,t) + t* Eij (coord) (i,j,t+1);
                }
            }
        }
    }
}

void integrator::setupEcd(){
    Ecd.set_size(3);
    int I,J,T;
    for(int coord=0;coord<3;coord++){
        I = pCijk(coord);
        J = pDijk(coord);
        T = pCijk(coord)+pDijk(coord);
        Ecd(coord).set_size(I+5,J+5,T+6); //not sure about these values (could be lower)
        Ecd(coord).zeros();
        //Ecd(coord) (1,1,1) = exp(-(c*d/q)*(Xcd(0)*Xcd(0)+Xcd(1)*Xcd(1)+Xcd(2)*Xcd(2)));
        Ecd(coord) (1,1,1) = exp(-(c*d/q)*(Xcd(coord)*Xcd(coord)));
        for(int i=1;i<I+4;i++){
            for(int t=1;t<T+5; t++){
                //Ecd(coord) (i+1,1,t) = Ecd(coord) (i,1,t-1)/(2*q) + Xqc(coord)*Ecd(coord) (i,1,t) + t* Eij (coord) (i,1,t+1);
                Ecd(coord) (i+1,1,t) = Ecd(coord) (i,1,t-1)/(2*q) + Xqc(coord)*Ecd(coord) (i,1,t) + t* Ecd (coord) (i,1,t+1); //maybe getting late ?
            }
        }
        for(int j=1;j<J+4;j++){
            for(int i=1;i<I+4;i++){
                for(int t=1;t<T+5; t++){
                    //Ecd(coord) (i,j+1,t) = Ecd(coord) (i,j,t-1)/(2*q) + Xqd(coord)*Ecd(coord) (i,j,t) + t* Eij (coord) (i,j,t+1);
                    Ecd(coord) (i,j+1,t) = Ecd(coord) (i,j,t-1)/(2*q) + Xqd(coord)*Ecd(coord) (i,j,t) + t* Ecd (coord) (i,j,t+1);
                }
            }
        }
    }
}

void integrator::setupRtuv(vec3 &nucleiPos){
    int T,U,V,N;
    T = pAijk(0)+pBijk(0);
    U = pAijk(1)+pBijk(1);
    V = pAijk(2)+pBijk(2);

    vec Aa = {T,U,V};

    N = T+U+V;
    //Rtuv = field<cube>(); //This is temporarily inserted to debug the integrator
    Rtuv.set_size(N+2);

    Rpc = P - nucleiPos;
    Rpc2 = Rpc(0)*Rpc(0)+Rpc(1)*Rpc(1)+Rpc(2)*Rpc(2);

    //BoysFunction boys(Aa.max()+1); //this function behaves oddly!!
    boys.setx(p*Rpc2);
    for(int n=0;n<N+2;n++){
        Rtuv(n).set_size(T+3,U+3,V+3);
        Rtuv(n).zeros();
        Rtuv(n) (1,1,1) = pow((-2.0*p),(double) n)*boys.returnValue(n);
    }
    for(int t=1;t<T+2;t++){
        for(int n=0;n<N+1;n++){
            Rtuv(n) (t+1,1,1) = (t-1)*Rtuv(n+1) (t-1,1,1) + Rpc(0) * Rtuv(n+1) (t,1,1);
        }
    }
    for(int u=1;u<U+2;u++){
        for(int t=1;t<T+2;t++){
            for(int n=0;n<N+1;n++){
                Rtuv(n) (t,u+1,1) = (u-1)*Rtuv(n+1) (t,u-1,1) + Rpc(1) * Rtuv(n+1) (t,u,1);
            }
        }
    }
    for(int v=1;v<V+2;v++){
        for(int u=1;u<U+2;u++){
            for(int t=1;t<T+2;t++){
                for(int n=0;n<N+1;n++){
                    Rtuv(n) (t,u,v+1) = (v-1)*Rtuv(n+1) (t,u,v-1) + Rpc(2) * Rtuv(n+1) (t,u,v);
                }
            }
        }
    }
}


void integrator::setupRtau(){
    //shared variables for A,B,C,D
    alpha = p*q/(p+q);
    Rpq = P-Q;

    int T,U,V,N,Tau,Ny,Phi;
    T = pAijk(0)+pBijk(0);
    U = pAijk(1)+pBijk(1);
    V = pAijk(2)+pBijk(2);
    Tau = pCijk(0)+pDijk(0);
    Ny  = pCijk(1)+pDijk(1);
    Phi = pCijk(2)+pDijk(2);

    //vec Aa = {T,U,V,Tau,Ny,Phi};

    N = T+U+V+Tau+Ny+Phi;
    //Rtau = field<cube>(); //This is temporarily inserted to debug the integrator
    Rtau.set_size(N+2);

    Rpq2 = Rpq(0)*Rpq(0)+Rpq(1)*Rpq(1)+Rpq(2)*Rpq(2);

    //BoysFunction boys(Aa.max()+1); //this function behaves oddly!!
    boys.setx(alpha*Rpq2);

    for(int n=0;n<N+2;n++){
        Rtau(n).set_size(T+Tau+3,U+Ny+3,V+Phi+3);
        Rtau(n).zeros();
        Rtau(n) (1,1,1) = pow((-2.0*alpha),(double) n)*boys.returnValue(n);
    }

    for(int t=1;t<T+Tau+2;t++){
        for(int n=0;n<N+1;n++){
            Rtau(n) (t+1,1,1) = (t-1)*Rtau(n+1) (t-1,1,1) + Rpq(0) * Rtau(n+1) (t,1,1);
        }
    }
    for(int u=1;u<U+Ny+2;u++){
        for(int t=1;t<T+Tau+2;t++){
            for(int n=0;n<N+1;n++){
                Rtau(n) (t,u+1,1) = (u-1)*Rtau(n+1) (t,u-1,1) + Rpq(1) * Rtau(n+1) (t,u,1);
            }
        }
    }
    for(int v=1;v<V+Phi+2;v++){
        for(int u=1;u<U+Ny+2;u++){
            for(int t=1;t<T+Tau+2;t++){
                for(int n=0;n<N+1;n++){
                    Rtau(n) (t,u,v+1) = (v-1)*Rtau(n+1) (t,u,v-1) + Rpq(2) * Rtau(n+1) (t,u,v);
                }
            }
        }
    }
}

double integrator::overlap(){
    //The overlap integral <pA|pB>
    double result = wA*wB*pow(sqrt(pi/p),3);
    int I,J;
    for(int coord=0; coord<3;coord++){
        I = pAijk(coord);
        J = pBijk(coord);
        result *= Eij(coord) (I+1,J+1,1);
    }
    return result;
}

double integrator::kinetic(){
    //The kinetic integral <pA|T|pB>
    S = sqrt(pi/p);
    for(int coord=0;coord<3;coord++){
        int i = pAijk(coord);
        int j = pBijk(coord);
        //double Sij_2 = Eij(coord) ((int) pAijk(coord)+1, (int) pBijk(coord)+1+2);
        Sijk(coord) = S*Eij(coord) ((int) pAijk(coord)+1, (int) pBijk(coord)+1, 1);
        Tijk(coord) = 4*b*b*S *Eij(coord) ((int) pAijk(coord)+1, (int) pBijk(coord)+3, 1) - 2*b*(2*pBijk(coord)+1)*S*Eij(coord) ((int) pAijk(coord)+1, (int) pBijk(coord)+1, 1);
        if(-1<pBijk(coord)-2){
            Tijk(coord) += pBijk(coord)*(pBijk(coord)-1)*S*Eij(coord) ((int) pAijk(coord)+1, (int) pBijk(coord)-1, 1);
            //cout << "AA" << endl;
        }
        //cout << Sijk(coord) << endl;
        //cout << Tijk(coord) << endl;
    }

    return -.5*wA*wB*(Tijk(0)*Sijk(1)*Sijk(2) + Tijk(1)*Sijk(2)*Sijk(0)+Tijk(2)*Sijk(0)*Sijk(1));
    //return -.5*wA*wB*(Tijk(0) + Tijk(1)+Tijk(2));
}

double integrator::pNuclei(){
    //The particle-nuclei interaction <pA|V|pB>
    int T,U,V;
    T = pAijk(0)+pBijk(0);
    U = pAijk(1)+pBijk(1);
    V = pAijk(2)+pBijk(2);
    double result = 0;
    for(int t=0;t<T+1;t++){
        for(int u=0;u<U+1;u++){
            for(int v=0;v<V+1;v++){
                result += Rtuv(0) (t+1,u+1,v+1)* Eij(0) ((int) pAijk(0)+1, (int) pBijk(0)+1, t+1)*Eij(1) ((int) pAijk(1)+1, (int) pBijk(1)+1, u+1)*Eij(2) ((int) pAijk(2)+1, (int) pBijk(2)+1, v+1);
            }
        }
    }
    return wA*wB*result*(2*pi/p);
}

double integrator::pp(Primitive &pC, Primitive &pD){
    //The particle-particle interaction <pA pB|r**-1 |pC pD>_AS
    c = pC.exponent();          // exponential constant.
    C = pC.nucleusPosition();   // nucleus pC position
    pCijk(0) = pC.xExponent();
    pCijk(1) = pC.yExponent();
    pCijk(2) = pC.zExponent();
    wC = pC.weight();

    d = pD.exponent();          // exponential constant.
    D = pD.nucleusPosition();   // nucleus pD position
    pDijk(0) = pD.xExponent();
    pDijk(1) = pD.yExponent();
    pDijk(2) = pD.zExponent();
    wD = pD.weight();

    //shared variables for C and D
    q = c+d;
    Q = (c*C+d*D)/q;
    Xcd = C-D;
    Xqc = Q-C;
    Xqd = Q-D;
    Xcd2 = Xcd(0)*Xcd(0)+Xcd(1)*Xcd(1)+Xcd(2)*Xcd(2);

    setupEcd();
    setupRtau();

    int T,U,V,Tau,Ny,Phi;
    T = pAijk(0)+pBijk(0);
    U = pAijk(1)+pBijk(1);
    V = pAijk(2)+pBijk(2);
    Tau = pCijk(0)+pDijk(0);
    Ny  = pCijk(1)+pDijk(1);
    Phi = pCijk(2)+pDijk(2);
    double E1,E2,R1,F2;
    double result = 0;
    for(int t=0;t<T+1;t++){
        for(int u=0;u<U+1;u++){
            for(int v=0;v<V+1;v++){
                for(int tau=0;tau<Tau+1;tau++){
                    for(int ny=0;ny<Ny+1;ny++){
                        for(int phi=0;phi<Phi+1;phi++){
                            E1 = Eij(0) ((int) pAijk(0)+1, (int) pBijk(0)+1, t+1)  *Eij(1) ((int) pAijk(1)+1, (int) pBijk(1)+1, u+1) *Eij(2) ((int) pAijk(2)+1, (int) pBijk(2)+1, v+1);
                            E2 = Ecd(0) ((int) pCijk(0)+1, (int) pDijk(0)+1, tau+1)*Ecd(1) ((int) pCijk(1)+1, (int) pDijk(1)+1, ny+1)*Ecd(2) ((int) pCijk(2)+1, (int) pDijk(2)+1, phi+1);
                            R1 = Rtau(0) (t+tau+1,u+ny+1,v+phi+1)*pow(-1,tau+ny+phi);
                            result += E1*E2*R1;
                        }
                    }
                }
            }
        }
    }
    return wA*wB*wC*wD*result*(2*pow(pi,2.5))/(p*q*sqrt(p+q));
}
