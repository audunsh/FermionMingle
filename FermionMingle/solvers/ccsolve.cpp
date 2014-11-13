#include "solvers/ccsolve.h"
#include "solvers/rhfsolve.h"
#include <iomanip>
#include "solvers/hfsolve.h"
ccsolve::ccsolve()
{
}

ccsolve::ccsolve(HFSolve object)
{
    /* --- CCSolve Main ---
     * This class should be initialized with a HFSolve object containing an energy-minimized basis for the Coupled-Cluster (CC) solver.
     * The CC Solver will then perform the following operations
     * 1. Set up a new basis from atomic to molecular orbitals
     * 2. Initialize the cluster amplitudes with a given initial value
     * 3. Solve the amplitude equations
     * 4. Return the CC-energy
    */
    hfobject = object;
    //nElectrons = nElect; //Fermi level is defined by number of electrons

    //retranslate();
    //CCSD();

    //cout << "Energy:" << energy(t2c, t1c) << endl;

}

void ccsolve::CombineUHFCoefficients(){
    mat C;
    C.zeros(nStates, 2*nStates);
    //C.zeros(2*nStates, 2*nStates);
    for(int i = 0; i<nStates; i++){
        //colu.subvec(0,nStates/2-1) = hfobject.Cu.col(i);

        //C.col(2*i).subvec(0,nStates/2) =  hfobject.Cu.col(i);
        //C.col(2*i + 1).subvec(nStates/2+1, nStates+1) = hfobject.Cd.col(i);

        C.col(2*i) =  hfobject.Cu.col(i);
        C.col(2*i + 1) = hfobject.Cd.col(i);
    }
    UHFC = C;
}

void ccsolve::init_UHF_basis(){
    nStates = hfobject.Cd.n_cols; // + hfobject.Cu.n_cols; //Double spin later
    CombineUHFCoefficients();
    UHFC.print();
    cout << " " << endl;
    SetupMinimizedUHFBasis();
    //vmin.print();
    //vminu.print();
    //Antisymmetrize();
    //ExpandMinimizedBasis();
    SetupT1();
    SetupT2();
    vmin.print();

}

void ccsolve::Antisymmetrize(){
    //nStates *= 2;
    temp_mo = vmin;
    vmin.set_size(nStates, nStates);
    fmin.set_size(nStates, nStates);
    fmin.zeros();
    //hfobject.C.col(4)*=-1;

    double val1 = 0.0;
    double val2 = 0.0;
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            //previously aibj
            fmin(p,q) = GetUncoupledElement(p,q);
            vmin(p,q) = zeros(nStates,nStates);
        }
    }
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            for(int r = 0; r<nStates; r++){
                for(int s=0; s<nStates; s++){

                    val1 = equalfunc(p%2,q%2) * equalfunc(r%2,s%2) * temp_mo(p/2,q/2)(r/2,s/2);
                    val2 = equalfunc(p%2,s%2) * equalfunc(r%2,q%2) * temp_mo(p/2,s/2)(r/2,q/2);

                    vmin(p,r)(q,s) = val1 - val2;



                    //vmin(p,r)(q,s) = hfobject.Bs.state(p,r,q,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,s/2)(r/2,q/2)); //THIS PRODUCES SAME RESULTS (vmin) AS FROM NORDLI
                    //vmin(p,r)(q,s) = hfobject.Bs.state(p,r,q,s, temp_mo(p/2,r/2)(q/2,s/2), temp_mo(p/2,r/2)(q/2,s/2)); //THIS PRODUCES SAME RESULTS (vmin) AS FROM NORDLI

                    //vmin(p,q)(r,s) = hfobject.Bs.state(p,q,r,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,q/2)(s/2,r/2)); //Trying to make minor changes in indexing
                    //vmin(p,q)(r,s) = hfobject.Bs.state(p,q,r,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,q/2)(s/2,r/2)); //Minor editing, producing errors
                }
            }
        }
    }
    cout << endl;
    //vmin.print();
    //vmin(13,12)(0,1)=0.0;
}

void ccsolve::init_RHF_basis(){
    nStates = hfobject.C.n_cols;
    //cout << "CCSolve Found " << nStates << " number of basis functions in the RHF Roothaan expansion." << endl;
    SetupMinimizedBasis();
    ExpandMinimizedBasis(); //Include spin orthogonality
    SetupT1();
    SetupT2();
    //vmin.print();

}

double ccsolve::GetUncoupledElement(int a, int b){
    double sm = 0;
    if(a==b){
        sm = hfobject.epsilon(a/2);
    }


    return sm;
}

double ccsolve::UHFCoefficient(int i,int j,int k,int l,int a,int b,int c,int d){

    //double val1 = equalfunc(a%2,b%2) * equalfunc(c%2,d%2) *equalfunc(i%2,j%2) * equalfunc(k%2,l%2) * hfobject.Bs.v(i/2,j/2)(k/2,l/2);
    //double val2 = equalfunc(a%2,d%2) * equalfunc(c%2,b%2) *equalfunc(i%2,l%2) * equalfunc(k%2,j%2) * hfobject.Bs.v(i/2,l/2)(k/2,j/2);

    //vmin(p,r)(q,s) = val1 - val2;

}



double ccsolve::GetCoupledElement(int a, int b, int c, int d){
    //Unsure about this transformation
    double sm = 0.0;
    for(int i=0; i<nStates; i++){
        for(int j=0; j<nStates; j++){
            for(int k=0; k<nStates; k++){
                for(int l=0; l<nStates; l++){
                    //sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.coupledMatrix(i,j)(k,l);
                    //sm += hfobject.C(a,i)*hfobject.C(b,j)*hfobject.C(c,k)*hfobject.C(d,l)*hfobject.Bs.v(i,j)(k,l);

                    sm += hfobject.C(i,a)*hfobject.C(j,b)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.Bs.v(i,j)(k,l);
                    //sm += hfobject.C(i,a)*hfobject.C(j,b)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.coupledMatrix(i,k)(j,l);
                    //sm += hfobject.C(i,a)*hfobject.C(j,b)*hfobject.C(k,c)*hfobject.C(l,d)*hfobject.coupledMatrix(i,j)(k,l);
                }
            }
        }
    }
    return sm;
}

double ccsolve::GetCoupledUHFElement(int a, int b, int c, int d){
    double sm = 0.0;
    for(int i=0; i<nStates/2; i++){
        for(int j=0; j<nStates/2; j++){
            for(int k=0; k<nStates/2; k++){
                for(int l=0; l<nStates/2; l++){
                     //sm += UHFC(i,a)*UHFC(j,b)*UHFC(k,c)*UHFC(l,d)*hfobject.Bs.v(i,j)(k,l);
                     sm += hfobject.Cu(i,a/2)*hfobject.Cu(j,b/2)*hfobject.Cu(k,c/2)*hfobject.Cu(l,d/2)*hfobject.Bs.v(i,j)(k,l);
                     sm += hfobject.Cd(i,a/2)*hfobject.Cd(j,b/2)*hfobject.Cd(k,c/2)*hfobject.Cd(l,d/2)*hfobject.Bs.v(i,j)(k,l);
                     sm +=4* hfobject.Cu(i,a/2)*hfobject.Cd(j,b/2)*hfobject.Cu(k,c/2)*hfobject.Cd(l,d/2)*hfobject.Bs.v(i,j)(k,l);

                }
            }
        }
    }
    return  sm; //equalfunc(a%2,c%2)*equalfunc(b%2,d%2)*
}

void ccsolve::SetupMinimizedUHFBasis(){
    //cout << "Setting up the minimized basis." << endl;
    //UHFC.print();
    nStates *= 2;
    vmin.set_size(nStates,nStates);
    fmin.set_size(nStates,nStates);
    fmin.zeros();
    for(int a=0; a<nStates; a++){
        for(int b=0; b<nStates; b++){
            //fmin(a,b) = GetUncoupledElement(a,b);
            vmin(a,b) = zeros<mat>(nStates,nStates);}}
    double di = 0;
    double ex = 0;
    for(int p=0; p<nStates; p++){
        //fmin(a,a) = hfobject.epsilon(a/2);
        for(int q=0; q<nStates; q++){
            for(int r=0; r<nStates; r++){
                for(int s=0; s<nStates; s++){
                    //vmin(a,b)(c,d) = hfobject.P(b,c)*hfobject.P(a,d)*hfobject.Bs.v(a,b)(c,d);
                    //di = equalfunc(p%2,q%2)*equalfunc(r%2,s%2)*GetCoupledUHFElement(p,q,r,s) ;
                    //ex = equalfunc(p%2,s%2)*equalfunc(r%2,q%2)*GetCoupledUHFElement(p,s,r,q) ;
                    //vmin(p,r)(q,s) =di-ex; // equalfunc(p%2,r%2)*equalfunc(q%2,s%2)*GetCoupledUHFElement(p,q,r,s) ; //Adjusting this parameter to compensate for notational differences
                    vmin(p,r)(q,s) = GetCoupledUHFElement(p,q,r,s);
                }
            }
        }
    }
    //vmin.print();
    //cout << "Finished setting up the minimized basis." << endl;
}

void ccsolve::SetupMinimizedUHFBasis2(){
    //cout << "Setting up the minimized basis." << endl;
    UHFC.print();
    vmin.set_size(nStates,nStates);
    vminu.set_size(nStates,nStates);
    fmin.set_size(nStates,nStates);
    fmin.zeros();
    for(int a=0; a<nStates; a++){
        for(int b=0; b<nStates; b++){
            //fmin(a,b) = GetUncoupledElement(a,b);
            vmin(a,b) = zeros<mat>(nStates,nStates);
            vminu(a,b) = zeros<mat>(nStates,nStates);
        }
    }
    for(int p=0; p<nStates; p++){
        //fmin(a,a) = hfobject.epsilon(a/2);
        for(int q=0; q<nStates; q++){
            for(int r=0; r<nStates; r++){
                for(int s=0; s<nStates; s++){
                    //vmin(a,b)(c,d) = hfobject.P(b,c)*hfobject.P(a,d)*hfobject.Bs.v(a,b)(c,d);
                    //vmin(p,q)(r,s) = GetCoupledUHFElement(p,q,r,s) ; //Adjusting this parameter to compensate for notational differences
                    vmin(p,q)(r,s) = GetCoupledUHFElement(p,q,r,s) ; //Adjusting this parameter to compensate for notational differences
                    vminu(p,q)(r,s) = GetCoupledUHFElement(p,q,r,s) ; //Adjusting this parameter to compensate for notational differences
                }
            }
        }
    }
    //cout << "Finished setting up the minimized basis." << endl;
}


void ccsolve::SetupMinimizedBasis(){
    //cout << "Setting up the minimized basis." << endl;
    //hfobject.C.print();
    vmin.set_size(nStates,nStates);
    fmin.set_size(nStates,nStates);
    fmin.zeros();
    for(int a=0; a<nStates; a++){
        for(int b=0; b<nStates; b++){
            //fmin(a,b) = GetUncoupledElement(a,b);
            vmin(a,b) = zeros<mat>(nStates,nStates);}}
    for(int a=0; a<nStates; a++){
        //fmin(a,a) = hfobject.epsilon(a/2);
        for(int b=0; b<nStates; b++){
            for(int c=0; c<nStates; c++){
                for(int d=0; d<nStates; d++){
                    //vmin(a,b)(c,d) = hfobject.P(b,c)*hfobject.P(a,d)*hfobject.Bs.v(a,b)(c,d);
                    vmin(a,b)(c,d) = GetCoupledElement(a,b,c,d) ; //Adjusting this parameter to compensate for notational differences
                }
            }
        }
    }
    //cout << "Finished setting up the minimized basis." << endl;
}

void ccsolve::ExpandMinimizedBasis(){
    nStates*= 2;
    temp_mo = vmin;
    vmin.set_size(nStates, nStates);
    fmin.set_size(nStates, nStates);
    fmin.zeros();
    //hfobject.C.col(4)*=-1;

    double val1 = 0.0;
    double val2 = 0.0;
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            //previously aibj
            fmin(p,q) = GetUncoupledElement(p,q);
            vmin(p,q) = zeros(nStates,nStates);
        }
    }
    for(int p = 0; p<nStates; p++){
        for(int q = 0; q<nStates; q++){
            for(int r = 0; r<nStates; r++){
                for(int s=0; s<nStates; s++){

                    val1 = equalfunc(p%2,q%2) * equalfunc(r%2,s%2) * temp_mo(p/2,q/2)(r/2,s/2);
                    val2 = equalfunc(p%2,s%2) * equalfunc(r%2,q%2) * temp_mo(p/2,s/2)(r/2,q/2);
                    vmin(p,r)(q,s) = val1 - val2;



                    //vmin(p,r)(q,s) = hfobject.Bs.state(p,r,q,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,s/2)(r/2,q/2)); //THIS PRODUCES SAME RESULTS (vmin) AS FROM NORDLI
                    //vmin(p,r)(q,s) = hfobject.Bs.state(p,r,q,s, temp_mo(p/2,r/2)(q/2,s/2), temp_mo(p/2,r/2)(q/2,s/2)); //THIS PRODUCES SAME RESULTS (vmin) AS FROM NORDLI

                    //vmin(p,q)(r,s) = hfobject.Bs.state(p,q,r,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,q/2)(s/2,r/2)); //Trying to make minor changes in indexing
                    //vmin(p,q)(r,s) = hfobject.Bs.state(p,q,r,s, temp_mo(p/2,q/2)(r/2,s/2), temp_mo(p/2,q/2)(s/2,r/2)); //Minor editing, producing errors
                }
            }
        }
    }
    cout << endl;
    //vmin.print();
    //vmin(13,12)(0,1)=0.0;
}

void ccsolve::expandC(){
    //expand the coefficient matrix to explicitly include spin dependence
    Cm.set_size(nStates,nStates);
    Cm.zeros();
    for(int i = 0; i<nStates/2; i++){
        for(int j = 0; j<nStates/2; j++){
            Cm(2*i,j) = hfobject.C(i,j);
            Cm(2*i+1,nStates/2+j) = hfobject.C(i,j);
        }
    }
    Cm.print();
}

double ccsolve::CCSD_Single(int a,int i){
    double S1 = 0;
    double S2a = 0;
    double S2b = 0;
    double S2c = 0;
    double S3a = 0;
    double S3b = 0;
    double S3c = 0;
    double S4a = 0;
    double S4b = 0;
    double S4c = 0;
    double S5a = 0;
    double S5b = 0;
    double S5c = 0;
    double S6 = 0;
    int k, l, c, d;

    S1 += fmin(a,i);

    for(k=0;k<nElectrons;k++){
        for(c=nElectrons;c<nStates;c++){
            S2a+=fmin(k,c)*t2(a,c)(i,k);
            for(d=nElectrons;d<nStates;d++){
                S2b+=vmin(a,k)(c,d)*t2c(c,d)(i,k);
            }
        }
    }


    for(k=0;k<nElectrons;k++){
        for(l=0;l<nElectrons;l++){
            for(c=nElectrons;c<nStates;c++){
                S2c-=vmin(k,l)(i,c)*t2c(a,c)(k,l);
            }
        }
    }

    for(c=nElectrons;c<nStates;c++){
        S3a+=fmin(a,c)*t1c(c,i);
    }


    for(k=0;k<nElectrons;k++){
        S3b-=fmin(k,i)*t1(a,k);
        for(c=nElectrons;c<nStates;c++){
            S3c += vmin(a,k)(i,c)*t1c(c,k);
        }
    }

    for(k=0;k<nElectrons;k++){
        for(l=0;l<nElectrons;l++){
            for(c=nElectrons;c<nStates;c++){
                S5c-=vmin(k,l)(i,c)*t1(a,k)*t1(c,l);
                for(d=nElectrons;d<nStates;d++){
                    S4a-=vmin(k,l)(c,d)*t1c(c,i)*t2c(a,d)(k,l);
                    S4b-=vmin(k,l)(c,d)*t1c(a,k)*t2c(c,d)(i,l);
                    S4c+=vmin(k,l)(c,d)*t1c(c,k)*t2c(d,a)(l,i);
                    S6-=vmin(k,l)(c,d)*t1c(c,i)*t1c(a,k)*t1c(d,l);
                }
            }
        }
    }

    for(k=0;k<nElectrons;k++){
        for(c=nElectrons;c<nStates;c++){
            S5a -= fmin(k,c)*t1c(c,i)*t1c(a,k);
            for(d=nElectrons; d<nStates; d++){
                S5b += vmin(a,k)(c,d)*t1c(c,i)*t1c(d,k);
            }
        }
    }

    //return S1+S2a+.5*S2b+.5*S2c+   0*S3a+0*S3b  +S3c+.5*S4a+.5*S4b+S4c+S5a+S5b+S5c+S6;
    return S1+S2a+.5*S2b+.5*S2c+   0*S3a+0*S3b  +S3c+.5*S4a+.5*S4b+S4c+S5a+S5b+S5c+S6; //HF-case
}

double ccsolve::CCSD_Double(int a, int b, int i, int j){
    int k,l,c,d;
    int nE = nElectrons;
    int nS = nStates;
    double D4a = 0;
    double D4b = 0;

    double D5a = 0;
    double D5b = 0;
    double D5c = 0;
    double D5d = 0;
    double D5e = 0;
    double D5f = 0;
    double D5g = 0;
    double D5h = 0;

    double D6a = 0;
    double D6b = 0;
    double D6c = 0;

    double D7a = 0;
    double D7b = 0;
    double D7c = 0;
    double D7d = 0;
    double D7e = 0;

    double D8a = 0;
    double D8b = 0;

    double D9  = 0;


    for(c=nE;c<nS;c++){
        D4a += vmin(a,b)(c,j)*t1c(c,i);
        D4a -= vmin(a,b)(c,i)*t1c(c,j);
    }

    for(k=0;k<nE;k++){
        D4b -= vmin(k,b)(i,j)*t1c(a,k);
        D4b += vmin(k,a)(i,j)*t1c(b,k);
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            D5a -= fmin(k,c)*t1c(c,i)*t2c(a,b)(k,j);
            D5a += fmin(k,c)*t1c(c,j)*t2c(a,b)(k,i);

            D5b -= fmin(k,c)*t1c(a,k)*t2c(c,b)(i,j);
            D5b += fmin(k,c)*t1c(b,k)*t2c(c,a)(i,j);

            for(d=nE;d<nS;d++){
                D5c += vmin(a,k)(c,d)*t1c(c,i)*t2c(d,b)(k,j);
                D5c -= vmin(a,k)(c,d)*t1c(c,j)*t2c(d,b)(k,i);
                D5c -= vmin(b,k)(c,d)*t1c(c,i)*t2c(d,a)(k,j);
                D5c += vmin(b,k)(c,d)*t1c(c,j)*t2c(d,a)(k,i);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                D5d -= vmin(k,l)(i,c)*t1c(a,k)*t2c(c,b)(l,j);
                D5d += vmin(k,l)(j,c)*t1c(a,k)*t2c(c,b)(l,i);
                D5d += vmin(k,l)(i,c)*t1c(b,k)*t2c(c,a)(l,j);
                D5d -= vmin(k,l)(j,c)*t1c(b,k)*t2c(c,a)(l,i);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D5e -= vmin(k,b)(c,d)*t1c(a,k)*t2c(c,d)(i,j);
                D5e += vmin(k,a)(c,d)*t1c(b,k)*t2c(c,d)(i,j);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                D5f += vmin(k,l)(c,j)*t1c(c,i)*t2c(a,b)(k,l);
                D5f -= vmin(k,l)(c,i)*t1c(c,j)*t2c(a,b)(k,l);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D5g += vmin(k,a)(c,d)*t1c(c,k)*t2c(d,b)(i,j);
                D5g -= vmin(k,b)(c,d)*t1c(c,k)*t2c(d,a)(i,j);
            }
        }
    }


    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                D5h -= vmin(k,l)(c,i)*t1c(c,k)*t2c(a,b)(l,j);
                D5h += vmin(k,l)(c,j)*t1c(c,k)*t2c(a,b)(l,i);
            }
        }
    }

    for(c=nE;c<nS;c++){
        for(d=nE;d<nS;d++){
            D6a += vmin(a,b)(c,d)*t1c(c,i)*t1c(d,j);
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            D6b += vmin(k,l)(i,j)*t1c(a,k)*t1c(b,l);
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            D6c -= vmin(k,b)(c,j)*t1c(c,i)*t1c(a,k);
            D6c += vmin(k,b)(c,i)*t1c(c,j)*t1c(a,k);
            D6c += vmin(k,a)(c,j)*t1c(c,i)*t1c(b,k);
            D6c -= vmin(k,a)(c,i)*t1c(c,j)*t1c(b,k);
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                for(d=nE; d<nS;d++){
                    D7a += vmin(k,l)(c,d)*t1c(c,i)*t1c(d,j)*t2c(a,b)(k,l);
                    D7b += vmin(k,l)(c,d)*t1c(a,k)*t1c(b,l)*t2c(c,d)(i,j);

                    D7c -= vmin(k,l)(c,d)*t1c(c,i)*t1c(a,k)*t2c(d,b)(l,j);
                    D7c += vmin(k,l)(c,d)*t1c(c,j)*t1c(a,k)*t2c(d,b)(l,i);
                    D7c += vmin(k,l)(c,d)*t1c(c,i)*t1c(b,k)*t2c(d,a)(l,j);
                    D7c -= vmin(k,l)(c,d)*t1c(c,j)*t1c(b,k)*t2c(d,a)(l,i);

                    D7d -= vmin(k,l)(c,d)*t1c(c,k)*t1c(d,i)*t2c(a,b)(l,j);
                    D7d += vmin(k,l)(c,d)*t1c(c,k)*t1c(d,j)*t2c(a,b)(l,i);

                    D7e -= vmin(k,l)(c,d)*t1c(c,k)*t1c(a,l)*t2c(d,b)(i,j);
                    D7e += vmin(k,l)(c,d)*t1c(c,k)*t1c(b,l)*t2c(d,a)(i,j);
                }
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D8a += vmin(k,b)(c,d)*t1c(c,i)*t1c(a,k)*t1c(d,j);
                D8a -= vmin(k,a)(c,d)*t1c(c,i)*t1c(b,k)*t1c(d,j);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(c=nE;c<nS;c++){
            for(d=nE;d<nS;d++){
                D8b += vmin(k,l)(c,j)*t1c(c,i)*t1c(a,k)*t1c(b,l);
                D8b -= vmin(k,l)(c,i)*t1c(c,j)*t1c(a,k)*t1c(b,l);
            }
        }
    }

    for(k=0;k<nE;k++){
        for(l=0;l<nE;l++){
            for(c=nE;c<nS;c++){
                for(d=nE; d<nS;d++){
                    D9 += vmin(k,l)(c,d)*t1c(c,i)*t1c(d,j)*t1c(a,k)*t1c(b,l);
                }
            }
        }
    }

    return D4a+D4b + D5a+D5b+D5c+D5d+.5*D5e+.5*D5f+ D5g+D5h + D6a + D6b + D6c + .5*D7a + .5*D7b + D7c + D7d + D7e + D8a + D8b + D9;
}

double ccsolve::CCDQ(int a, int b, int i, int j){
    //Quadratic contributions to the CCD amplitude equation
    double Qa = 0;
    double Qb = 0;
    double Qc = 0;
    double Qd = 0;
    //int a,b,i,j; //These should be included in the function call
    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qa += vmin(k,l)(c,d)*t2c(c,d)(i,j)*t2c(a,b)(k,l);
                }
            }
        }
    }
    //Qa/=4.0;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    //Qb += vmin(k,l)(c,d)*t2(a,c)(i,k)*t2(d,b)(l,j); (*) Conflictin notation in SB
                    //Qb -= vmin(k,l)(c,d)*t2(a,c)(j,k)*t2(d,b)(l,i);

                    Qb += vmin(k,l)(c,d)*t2c(a,c)(i,k)*t2c(b,d)(j,l);
                    Qb -= vmin(k,l)(c,d)*t2c(a,c)(j,k)*t2c(b,d)(i,l);
                }
            }
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qc -= vmin(k,l)(c,d)*t2c(c,d)(k,i)*t2c(a,b)(l,j);
                    Qc += vmin(k,l)(c,d)*t2c(c,d)(k,j)*t2c(a,b)(l,i);
                }
            }
        }
    }
    //Qc *= .5;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qd -= vmin(k,l)(c,d)*t2c(c,a)(k,l)*t2c(d,b)(i,j);
                    Qd += vmin(k,l)(c,d)*t2c(c,b)(k,l)*t2c(d,a)(i,j);
                }
            }
        }
    }
    //Qd *= .5;
    return 0.25*Qa + Qb + 0.5*(Qc + Qd);
}

double ccsolve::CCDQ2(int a, int b, int i, int j, field<mat> tf, double Qac, double Qbc, double Qcc, double Qdc){
    //Quadratic contributions to the CCD amplitude equation
    //Call the elements in the following order: vmin(i,j)(a,b)
    double Qa = 0;
    double Qb = 0;
    double Qc = 0;
    double Qd = 0;
    //int a,b,i,j; //These should be included in the function call
    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qa += vmin(k,l)(c,d)*tf(c,d)(i,j)*tf(a,b)(k,l);
                }
            }
        }
    }
    //Qa/=4.0;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    //Qb += vmin(k,l)(c,d)*t2(a,c)(i,k)*t2(d,b)(l,j); (*) Conflictin notation in SB
                    //Qb -= vmin(k,l)(c,d)*t2(a,c)(j,k)*t2(d,b)(l,i);

                    Qb += vmin(k,l)(c,d)*tf(a,c)(i,k)*tf(b,d)(j,l);
                    Qb -= vmin(k,l)(c,d)*tf(a,c)(j,k)*tf(b,d)(i,l);
                }
            }
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qc -= vmin(k,l)(c,d)*tf(d,c)(i,k)*tf(a,b)(l,j);
                    Qc += vmin(k,l)(c,d)*tf(d,c)(j,k)*tf(a,b)(l,i);
                }
            }
        }
    }
    //Qc *= .5;

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            for(int c = nElectrons; c<nStates; c++){
                for(int d = nElectrons; d<nStates; d++){
                    Qd -= vmin(k,l)(c,d)*tf(a,c)(l,k)*tf(d,b)(i,j);
                    Qd += vmin(k,l)(c,d)*tf(b,c)(l,k)*tf(d,a)(i,j);
                }
            }
        }
    }
    //Qd *= .5;
    //return 0.25*Qa + Qb + 0.5*(Qc + Qd);
    return Qac*Qa + Qbc*Qb + Qcc*Qc + Qdc*Qd;
}

double ccsolve::CCDL2(int a, int b, int i, int j, field<mat> tf, double L1ac, double L1bc,double L2ac,double L2bc,double L2cc){
    //Linear contributions to the CCD amplitude equation.
    //Call the elements in the following order: v(i,j)(a,b)
    double L1a = 0.0;
    double L1b = 0.0;
    double L2a = 0.0;
    double L2b = 0.0;
    double L2c = 0.0;

/*
    for(int c = nElectrons; c<nStates; c++){
        L1a += fmin(b,c)*tf(a,c)(i,j) - fmin(a,c)*tf(b,c)(i,j);
    }

    for(int k = 0; k<nElectrons; k++){
        L1b -= fmin(k,j)*tf(a,b)(i,k);
        L1b += fmin(k,i)*tf(a,b)(j,k);
    }
    */

    for(int c = nElectrons; c<nStates; c++){
        for(int d = nElectrons; d<nStates; d++){
            L2a += vmin(a,b)(c,d)*tf(c,d)(i,j); //H2O is highly sensitive to this contribution (2/10 2014)

            //L2a += vmin(a,c)(b,d)*t2(c,d)(i,j); //Shuffling some indices
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            //L2b += vmin(i,j)(k,l)*t2(a,b)(k,l); //(*) Conflicting notation in SB
            L2b += vmin(k,l)(i,j)*tf(a,b)(k,l);

            //L2b += vmin(k,i)(l,j)*t2(a,b)(k,l); //shuffling some indices
        }
    }

    for(int k= 0; k<nElectrons; k++){
        for(int c = nElectrons; c<nStates; c++){
            L2c += vmin(k,b)(c,j)*tf(a,c)(i,k);
            L2c -= vmin(k,a)(c,j)*tf(b,c)(i,k);
            L2c -= vmin(k,b)(c,i)*tf(a,c)(j,k);
            L2c += vmin(k,a)(c,i)*tf(b,c)(j,k);

            //L2c += vmin(b,k)(c,j)*t2c(a,c)(i,k);
            //L2c -= vmin(b,k)(c,i)*t2c(a,c)(j,k);
            //L2c -= vmin(a,k)(c,j)*t2c(b,c)(i,k);
            //L2c += vmin(a,k)(c,i)*t2c(b,c)(j,k);
        }
    }
    return L1ac*L1a + L1bc*L1b + L2ac*L2a + L2bc*L2b + L2cc*L2c; //0*.65*L2c;

}

double ccsolve::CCD_comp(int a, int b, int i, int j, field <mat> tf){
    double CC0_0a = 0.0;
    for(int c = nElectrons; c < nStates; c ++){
        CC0_0a += fmin(a,c)*tf(c,b)(i,j)-(fmin(b,c)*tf(c,a)(i,j));
    }
    CC0_0a *= -1.000000;


    double CC1_0a = 0.0;
    for(int k = 0; k < nElectrons; k ++){
        CC1_0a += fmin(k,j)*tf(a,b)(i,k)-(fmin(k,i)*tf(a,b)(j,k));
    }
    CC1_0a *= 1.000000;


    double CC4_0a = 0.0;
    for(int d = nElectrons; d < nStates; d ++){
        for(int c = nElectrons; c < nStates; c ++){
            CC4_0a += vmin(a,b)(c,d)*tf(c,d)(i,j);
        }
    }
    CC4_0a *= 0.500000;


    double CC5_0a = 0.0;
    for(int l = 0; l < nElectrons; l ++){
        for(int k = 0; k < nElectrons; k ++){
            CC5_0a += vmin(k,l)(i,j)*tf(a,b)(k,l);
        }
    }
    CC5_0a *= 0.500000;


    double CC6_0a = 0.0;
    for(int k = 0; k < nElectrons; k ++){
        for(int c = nElectrons; c < nStates; c ++){
            CC6_0a += vmin(a,k)(c,j)*tf(c,b)(i,k)-(vmin(b,k)(c,j)*tf(c,a)(i,k))-(vmin(a,k)(c,i)*tf(c,b)(j,k)-(vmin(b,k)(c,i)*tf(c,a)(j,k)));
        }
    }
    CC6_0a *= -1.000000;


    double CC12_1a = 0.0;
    for(int l = 0; l < nElectrons; l ++){
        for(int k = 0; k < nElectrons; k ++){
            for(int d = nElectrons; d < nStates; d ++){
                for(int c = nElectrons; c < nStates; c ++){
                    CC12_1a += vmin(k,l)(c,d)*tf(c,d)(i,k)*tf(a,b)(j,l)-(vmin(k,l)(c,d)*tf(c,d)(j,k)*tf(a,b)(i,l));
                }
            }
        }
    }
    CC12_1a *= -0.500000;


    double CC12_1b = 0.0;
    for(int l = 0; l < nElectrons; l ++){
        for(int k = 0; k < nElectrons; k ++){
            for(int d = nElectrons; d < nStates; d ++){
                for(int c = nElectrons; c < nStates; c ++){
                    CC12_1b += vmin(k,l)(c,d)*tf(c,d)(i,j)*tf(a,b)(k,l);
                }
            }
        }
    }
    CC12_1b *= 0.250000;


    double CC12_1c = 0.0;
    for(int d = nElectrons; d < nStates; d ++){
        for(int l = 0; l < nElectrons; l ++){
            for(int k = 0; k < nElectrons; k ++){
                for(int c = nElectrons; c < nStates; c ++){
                    CC12_1c += vmin(k,l)(c,d)*tf(c,a)(k,l)*tf(d,b)(i,j)-(vmin(k,l)(c,d)*tf(c,b)(k,l)*tf(d,a)(i,j));
                }
            }
        }
    }
    CC12_1c *= -0.500000;


    double CC12_1d = 0.0;
    for(int l = 0; l < nElectrons; l ++){
        for(int d = nElectrons; d < nStates; d ++){
            for(int k = 0; k < nElectrons; k ++){
                for(int c = nElectrons; c < nStates; c ++){
                    CC12_1d += vmin(k,l)(c,d)*tf(c,a)(i,k)*tf(d,b)(j,l)-(vmin(k,l)(c,d)*tf(c,b)(i,k)*tf(d,a)(j,l))-(vmin(k,l)(c,d)*tf(c,a)(j,k)*tf(d,b)(i,l)-(vmin(k,l)(c,d)*tf(c,b)(j,k)*tf(d,a)(i,l)));
                }
            }
        }
    }
    CC12_1d *= 0.500000;

    return CC12_1d + CC12_1c + CC12_1b + CC12_1a + CC6_0a + CC5_0a + CC4_0a + 0*CC1_0a + 0*CC0_0a;
}

double ccsolve::CCDL(int a, int b, int i, int j){
    //Linear contributions to the CCD amplitude equation.
    double L1a = 0.0;
    double L1b = 0.0;
    double L2a = 0.0;
    double L2b = 0.0;
    double L2c = 0.0;


    for(int c = nElectrons; c<nStates; c++){
        L1a += fmin(b,c)*t2c(a,c)(i,j) - fmin(a,c)*t2c(b,c)(i,j);
    }

    for(int k = 0; k<nElectrons; k++){
        L1b -= fmin(k,j)*t2c(a,b)(i,k);
        L1b += fmin(k,i)*t2c(a,b)(j,k);
    }

    for(int c = nElectrons; c<nStates; c++){
        for(int d = nElectrons; d<nStates; d++){
            L2a += vmin(a,b)(c,d)*t2c(c,d)(i,j); //H2O is highly sensitive to this contribution (2/10 2014)

            //L2a += vmin(a,c)(b,d)*t2(c,d)(i,j); //Shuffling some indices
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int l = 0; l<nElectrons; l++){
            //L2b += vmin(i,j)(k,l)*t2(a,b)(k,l); //(*) Conflicting notation in SB
            L2b += vmin(k,l)(i,j)*t2c(a,b)(k,l);

            //L2b += vmin(k,i)(l,j)*t2(a,b)(k,l); //shuffling some indices
        }
    }

    for(int k = 0; k<nElectrons; k++){
        for(int c = nElectrons; c<nStates; c++){
            L2c -= vmin(a,k)(c,j)*t2c(c,b)(i,k);
            L2c += vmin(a,k)(c,i)*t2c(c,b)(j,k);
            L2c += vmin(b,k)(c,j)*t2c(c,a)(i,k);
            L2c -= vmin(b,k)(c,i)*t2c(c,a)(j,k);

            //L2c += vmin(k,b)(c,j)*t2(a,c)(i,k);
            //L2c -= vmin(k,b)(c,i)*t2(a,c)(j,k);
            //L2c -= vmin(k,a)(c,j)*t2(b,c)(i,k);
            //L2c += vmin(k,a)(c,i)*t2(b,c)(j,k);

        }
    }
    //return 0*(L1a + L1b) + .5*(L2a + L2b)  + L2c; //0*.65*L2c;
    return 0*L1a + 0*L1b + .5*L2a + .5*L2b  + L2c; //0*.65*L2c;

}

void ccsolve::CCSD(int N_electrons){
    nElectrons = N_electrons;
    eprev = 100.0;
    maxiter = 100;
    counter = 0;
    //Set up and solve the CCSD-equation

    //cout << "Performing CCSD calculation." << endl;
    initT1(); //initial guess for the amplitudes
    initT2(); //Set up initial guess following S-B, p.289
    //t20 = t2new;
    double CCLinear, CCQuadratic;
    //cout << "MBPT2:" << energy(t2c,t1c) << endl;

    while(unconverged(.00000001)){
        //cout <<"::Current energy: " << energy(t2new,t1new) << endl;
        eprev = energy(t2c,t1c); //updating the energy
        //ReportEnergy();

        for(int b = nElectrons; b<nStates; b++){
            for(int a = b+1; a<nStates; a++){
                for(int j = 0; j<nElectrons; j++){
                    for(int i=j+1; i<nElectrons; i++){
                        //CCLinear = CCDL2(a,b,i,j, t2c,0.0,0.0,0.5,  0.5,  1.0);
                        //CCQuadratic = CCDQ2(a,b,i,j, t2c,0.25,1.0,0.5,0.5);

                        CCLinear = CCD_comp(a,b,i,j, t2c); //Experimenting with code generation from CCAlgebra, proof-of-consept: CCD generated equations working!
                        CCQuadratic = 0;


                        //t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCSD_Double(a,b,i,j) +  CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Working
                        t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCSD_Double(a,b,i,j) +  CCLinear + CCQuadratic)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Refining
                        t2new(b,a)(i,j) = -t2new(a,b)(i,j);
                        t2new(a,b)(j,i) = -t2new(a,b)(i,j);
                        t2new(b,a)(j,i) =  t2new(a,b)(i,j);

                    }
                }
            }
        }
        for(int a=nElectrons;a<nStates;a++){
            for(int i=0;i<nElectrons;i++){
                t1new(a,i) = CCSD_Single(a,i)/(fmin(i,i)-fmin(a,a));
            }
        }
        t2c = t2new;
        t1c = t1new;

    }
    correlation_energy = energy(t2c,t1c);
}

void ccsolve::retranslate(){
    temp_mo = vmin;
    int a,b,c,d;
    for(a=0;a<nStates;a++){
        for(b=0;b<nStates;b++){
            for(c=0;c<nStates;c++){
                for(d=0;d<nStates;d++){
                    vmin(a,b)(c,d) = temp_mo(a,d)(c,b);
                }
            }
        }
    }
}

void ccsolve::ReportEnergy(){
    int naW = 14;
    if(counter == 1){
        cout << "----------------------------------" << endl;
        cout << "Iteration     Correlation" << endl;
        cout << "---------------------------------" << endl;
    }
    cout << left << setw(naW) << setfill(' ') << counter;
    cout << left << setw(naW) << setfill(' ') << eprev ; //printing out results
    cout << endl;
}

void ccsolve::CCD(){
    //retranslate();
    eprev = 1.0;
    maxiter = 100;
    counter = 0;
    double CCDL2c = 0.0; //contributions
    double CCDQ2c = 0.0;


    //Set up and solve the CCSD-equation
    cout << "Performing CCD calculation." << endl;
    initT2(); //Set up initial guess following S-B, p.289


    cout << endl;
    //vmin.print();

    //cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    //cout << "Iteration     CCDEnergy     MBPT(2)        Linear        L2a           L2b           L2c           Quadratic     Qa            Qb            Qc            Qd" << endl;
    //cout << "-------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    while(unconverged(.00000001)){
        CCDL2c = 0.0; //contributions
        CCDQ2c = 0.0;

        eprev = energy(t2c, t1c); //updating the energy
        ReportEnergy();
        //cout << left << setw(naW) << setfill(' ') << counter;
        //cout << left << setw(naW) << setfill(' ') << eprev ; //printing out results

        //cout << left << setw(naW) << setfill(' ') << energy(t2new0, t1c);
        //cout << left << setw(naW) << setfill(' ') << energy(t2newL, t1c); //printing out results
        //cout << left << setw(naW) << setfill(' ') << energy(t2newL2a, t1c); //printing out results
        //cout << left << setw(naW) << setfill(' ') << energy(t2newL2b, t1c); //printing out results
        //cout << left << setw(naW) << setfill(' ') << energy(t2newL2c, t1c); //printing out results

        //cout << left << setw(naW) << setfill(' ') << energy(t2newQ, t1c); //printing out results
        //cout << left << setw(naW) << setfill(' ') << energy(t2newQa, t1c); //printing out results
        //cout << left << setw(naW) << setfill(' ') << energy(t2newQb, t1c); //printing out results
        //cout << left << setw(naW) << setfill(' ') << energy(t2newQc, t1c); //printing out results
        //cout << left << setw(naW) << setfill(' ') << energy(t2newQd, t1c); //printing out results



        for(int b = nElectrons; b<nStates; b++){
            for(int a = b+1; a<nStates; a++){
                for(int j = 0; j<nElectrons; j++){
                    for(int i=j+1; i<nElectrons; i++){

                        //t2new(a,b)(i,j) = (vmin(i,j)(a,b) + CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Original, working until 1. iteration
                        //t2new(a,b)(i,j) = (vmin(i,j)(a,b) + CCDL2(a,b,i,j) + CCDQ2(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        //t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCDL2(a,b,i,j) + 0*CCDQ2(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));

                        t2new(a,b)(i,j) = (vmin(a,b)(i,j) + CCDL2(a,b,i,j, t2c,0.0,0.0,0.5,  0.5,  1.0) + CCDQ2(a,b,i,j, t2p,0.25,1.0,0.5,0.5))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2new(b,a)(i,j) = -t2new(a,b)(i,j);
                        t2new(a,b)(j,i) = -t2new(a,b)(i,j);
                        t2new(b,a)(j,i) =  t2new(a,b)(i,j);

                        t2new0(a,b)(i,j) = vmin(a,b)(i,j)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2new0(b,a)(i,j) = -t2new0(a,b)(i,j);
                        t2new0(a,b)(j,i) = -t2new0(a,b)(i,j);
                        t2new0(b,a)(j,i) =  t2new0(a,b)(i,j);

                        t2newL(a,b)(i,j) = CCDL2(a,b,i,j, t2c,0.0,0.0,0.5,0.5,1.0)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newL(b,a)(i,j) = -t2newL(a,b)(i,j);
                        t2newL(a,b)(j,i) = -t2newL(a,b)(i,j);
                        t2newL(b,a)(j,i) =  t2newL(a,b)(i,j);

                        t2newL2a(a,b)(i,j) = CCDL2(a,b,i,j, t2newL2ap,0.0,0.0,0.5,0,0)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newL2a(b,a)(i,j) = -t2newL2a(a,b)(i,j);
                        t2newL2a(a,b)(j,i) = -t2newL2a(a,b)(i,j);
                        t2newL2a(b,a)(j,i) =  t2newL2a(a,b)(i,j);

                        t2newL2b(a,b)(i,j) = CCDL2(a,b,i,j, t2newL2bp,0.0,0.0,0,0.5,0)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newL2b(b,a)(i,j) = -t2newL2b(a,b)(i,j);
                        t2newL2b(a,b)(j,i) = -t2newL2b(a,b)(i,j);
                        t2newL2b(b,a)(j,i) =  t2newL2b(a,b)(i,j);

                        t2newL2c(a,b)(i,j) = CCDL2(a,b,i,j, t2newL2cp,0.0,0.0,0,0,1.0)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newL2c(b,a)(i,j) = -t2newL2c(a,b)(i,j);
                        t2newL2c(a,b)(j,i) = -t2newL2c(a,b)(i,j);
                        t2newL2c(b,a)(j,i) =  t2newL2c(a,b)(i,j);


                        t2newQ(a,b)(i,j) = CCDQ2(a,b,i,j, t2newQp,0.25,1.0,0.5,0.5)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newQ(b,a)(i,j) = -t2newQ(a,b)(i,j);
                        t2newQ(a,b)(j,i) = -t2newQ(a,b)(i,j);
                        t2newQ(b,a)(j,i) =  t2newQ(a,b)(i,j);

                        t2newQa(a,b)(i,j) = CCDQ2(a,b,i,j, t2newQap,0.25,0,0,0)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newQa(b,a)(i,j) = -t2newQa(a,b)(i,j);
                        t2newQa(a,b)(j,i) = -t2newQa(a,b)(i,j);
                        t2newQa(b,a)(j,i) =  t2newQa(a,b)(i,j);

                        t2newQb(a,b)(i,j) = CCDQ2(a,b,i,j, t2newQbp,0,1.0,0,0)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newQb(b,a)(i,j) = -t2newQb(a,b)(i,j);
                        t2newQb(a,b)(j,i) = -t2newQb(a,b)(i,j);
                        t2newQb(b,a)(j,i) =  t2newQb(a,b)(i,j);

                        t2newQc(a,b)(i,j) = CCDQ2(a,b,i,j, t2newQcp,0,0,0.5,0)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newQc(b,a)(i,j) = -t2newQc(a,b)(i,j);
                        t2newQc(a,b)(j,i) = -t2newQc(a,b)(i,j);
                        t2newQc(b,a)(j,i) =  t2newQc(a,b)(i,j);

                        t2newQd(a,b)(i,j) = CCDQ2(a,b,i,j, t2newQdp,0,0,0,0.5)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));
                        t2newQd(b,a)(i,j) = -t2newQd(a,b)(i,j);
                        t2newQd(a,b)(j,i) = -t2newQd(a,b)(i,j);
                        t2newQd(b,a)(j,i) =  t2newQd(a,b)(i,j);




                        //t2new(i,j)(a,b) = (vmin(a,b)(i,j) + CCDL(a,b,i,j) + CCDQ(a,b,i,j))/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Original, working until 1. iteration
                        //t2new(i,j)(b,a) = -t2new(i,j)(a,b);
                        //t2new(j,i)(a,b) = -t2new(i,j)(a,b);
                        //t2new(j,i)(b,a) =  t2new(i,j)(a,b);
                    }
                }
            }
        }


        cout << endl;

        //updating the amplitudes
        //t2p = t2c;
        t2c = t2new;

        t2newQp = t2newQ;
        t2newLp = t2newL;

        t2newQap = t2newQa;
        t2newQbp = t2newQb;
        t2newQcp = t2newQc;
        t2newQdp = t2newQd;

        t2newL2ap = t2newL2a;
        t2newL2bp = t2newL2b;
        t2newL2cp = t2newL2c;


    }
}


double ccsolve::energy(field<mat> tf, mat t1f){
    double dE1 = 0;
    double dE2 = 0;
    double dE3 = 0;
    for(int i = 0; i<nElectrons; i++){
        for(int j = 0; j<nElectrons; j++){
            for(int a=nElectrons; a<nStates; a++){
                for(int b=nElectrons; b<nStates; b++){
                    //dE += .25*vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    //dE +=  .5*vmin(i,j)(a,b)*t1new(a,i)*t1new(b,j); //originals, working

                    //dE += .25*vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    //dE +=  .5*vmin(i,j)(a,b)*t1new(a,i)*t1new(b,j);

                    dE1 += vmin(i,j)(a,b)*tf(a,b)(i,j);
                    dE2 += vmin(i,j)(a,b)*t1f(a,i)*t1f(b,j);
                }
            }
        }
    }

    for(int i = 0; i<nElectrons; i++){
        for(int a=nElectrons; a<nStates; a++){
            //since i!=a this will never contribute to the energy when using a HF SD.
            dE3 += fmin(i,a)*t1f(a,i);
        }
    }

    return .25*dE1 + .5*dE2 + dE3;
}



double ccsolve::equalfunc(int a, int b){
    double ret = 0.0;
    if(a==b){
        ret = 1.0;
    }
    return ret;
}


double ccsolve::CCDenergy(field<mat>tf){
    double dE = 0;
    for(int i = 0; i<nElectrons; i++){
        for(int j = 0; j<nElectrons; j++){
            for(int a=nElectrons; a<nStates; a++){
                for(int b=nElectrons; b<nStates; b++){
                    //dE += vmin(i,j)(a,b)*t2new(a,b)(i,j);
                    dE += vmin(i,j)(a,b)*tf(a,b)(i,j);
                    //dE += .25*vmin(a,b)(i,j)*t2c(a,b)(i,j);


                }
            }
        }
    }
    return .25*dE;
}


bool ccsolve::unconverged(double tolerance){
    //Convergence criteria for CC iterations
    bool condition = true;
    double ec = energy(t2c,t1c);
    if(abs(ec-eprev)<tolerance){
        //cout << "CCsolve:Converged by surpassing tolerance in " << counter << " iterations." << endl;
        condition = false;
    }
    if(ec!=ec){condition = false;}
    if(counter > maxiter){
        //cout << "Maximum number of iterations exceeded. (" << maxiter << ")" << endl;
        condition = false;}
    counter += 1;
    return condition;
}


void ccsolve::initT2(){
    //initial values for T2 amplitudes
    for(int a=nElectrons; a<nStates; a++){
        for(int b=nElectrons; b<nStates; b++){
            for(int i=0;i<nElectrons;i++){
                for(int j=0;j<nElectrons;j++){
                    //t2new(a,b)(i,j) = vmin(i,j)(a,b)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b)); //Original, working
                    t2c(a,b)(i,j) = vmin(a,b)(i,j)/(fmin(i,i) + fmin(j,j) - fmin(a,a) - fmin(b,b));



                }
            }
        }
    }
    t2new = t2c;
    t2newLp   = t2c;
    t2newL2ap = t2c;
    t2newL2bp = t2c;
    t2newL2cp = t2c;

    t2newQp  = t2c;
    t2newQap = t2c;
    t2newQbp = t2c;
    t2newQcp = t2c;
    t2newQdp = t2c;


}


void ccsolve::initT1(){
    //initial value for T1 amplitudes
    for(int a=nElectrons; a<nStates; a++){
        for(int i=0; i<nElectrons; i++){
            t1new(i,a) = 0;//fmin(a,i)/(fmin(i,i) - fmin(a,a));
            t1c(i,a) = 0;//fmin(a,i)/(fmin(i,i) - fmin(a,a));
        }
    }
}


void ccsolve::SetupT1(){
    t1.set_size(nStates, nStates);
    t1.zeros();
    t1c.set_size(nStates, nStates);
    t1c.zeros();
    t1p.set_size(nStates, nStates);
    t1p.zeros();
    t1new.set_size(nStates, nStates);
    t1new.zeros();

    //cout << "Initialized t1 amplitudes." << endl;
}



void ccsolve::SetupT2(){
    t2.set_size(nStates, nStates);
    t2c.set_size(nStates, nStates);
    t2p.set_size(nStates, nStates);
    t2new.set_size(nStates, nStates);
    t2new0.set_size(nStates, nStates);

    t2newL.set_size(nStates, nStates);
    t2newL2a.set_size(nStates, nStates);
    t2newL2b.set_size(nStates, nStates);
    t2newL2c.set_size(nStates, nStates);

    t2newQ.set_size(nStates, nStates);
    t2newQa.set_size(nStates, nStates);
    t2newQb.set_size(nStates, nStates);
    t2newQc.set_size(nStates, nStates);
    t2newQd.set_size(nStates, nStates);

    for(int i = 0; i<nStates; i++){
        for(int j =0;j<nStates; j++){
            t2(i,j).set_size(nStates, nStates);
            t2(i,j).zeros();
            t2new(i,j).set_size(nStates, nStates);
            t2new(i,j).zeros();

            t2new0(i,j).set_size(nStates, nStates);
            t2new0(i,j).zeros();

            t2newL(i,j).set_size(nStates, nStates);
            t2newL(i,j).zeros();

            t2newL2a(i,j).set_size(nStates, nStates);
            t2newL2a(i,j).zeros();
            t2newL2b(i,j).set_size(nStates, nStates);
            t2newL2b(i,j).zeros();
            t2newL2c(i,j).set_size(nStates, nStates);
            t2newL2c(i,j).zeros();

            t2newQ(i,j).set_size(nStates, nStates);
            t2newQ(i,j).zeros();

            t2newQa(i,j).set_size(nStates, nStates);
            t2newQa(i,j).zeros();
            t2newQb(i,j).set_size(nStates, nStates);
            t2newQb(i,j).zeros();
            t2newQc(i,j).set_size(nStates, nStates);
            t2newQc(i,j).zeros();
            t2newQd(i,j).set_size(nStates, nStates);
            t2newQd(i,j).zeros();

            t2c(i,j).set_size(nStates, nStates);
            t2c(i,j).zeros();
            t2p(i,j).set_size(nStates, nStates);
            t2p(i,j).zeros();

        }
    }
    //cout << "Initialized t2 amplitudes." << endl;
}




