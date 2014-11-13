#include "solvers/uhfsolve.h"
#include "basis/basis.h"

uhfsolve::uhfsolve(){}

uhfsolve::uhfsolve(basis BS, int N_spin_up, int N_spin_down){
    //initialize solver with given basis, number of electrons (and now superfluos number of protons)
    //This solver follows the algo described in Thijssen, p74-76
    Bs = BS;
    nStates = Bs.Nstates;        //set number of states in basis

    nElectrons = N_spin_up + N_spin_down;              //set number of electrons
    nElectronsU = N_spin_up;
    nElectronsD = N_spin_down;
    //nProtons = Z;                //set number of protons, may be removed

    //initializing all matrices and vectors
    C.zeros(nStates,nElectrons/2);//set initial C equal to the unit matrix

    //Cu.zeros(nStates,nElectronsU);//set initial C equal to the unit matrix
    //Cd.zeros(nStates,nElectronsD);//set initial C equal to the unit matrix

    Cu.zeros(nStates,nStates);//set initial C equal to the unit matrix
    Cd.zeros(nStates,nStates);//set initial C equal to the unit matrix
    //Cu.randn(nStates, nStates);
    //Cd.randn(nStates, nStates);


    F.zeros(nStates,nStates);     //initialize Fock matrix
    Fu.zeros(nStates,nStates);     //initialize Fock matrix
    Fd.zeros(nStates,nStates);     //initialize Fock matrix

    P.zeros(nStates,nStates);     //initialize Density matrix
    Pu.zeros(nStates,nStates);     //initialize Density matrix
    Pd.zeros(nStates,nStates);     //initialize Density matrix

    U.zeros(nStates,nStates);     //initialize Unitary matrix
    G.zeros(nStates,nStates);     //initialize Fock-component matrix
    Fprime.zeros(nStates,nStates);//transformed Fock matrix
    epsilon.zeros(nStates);       //eigenvalues from diagonalization
    epsilon_prev.zeros(nStates);  //eigenvalues from previous diagonalization

    setupCoupledMatrix();  //import particle-particle interaction integrals
    //setupCoupledMatrix_unused();  //import particle-particle interaction integrals
    setupP();              //set up initial density matrix
    s_diag.zeros(nStates);
}

void uhfsolve::reset(basis BS, int N, int Z){
    //initialize solver with given basis, number of electrons (and now superfluos number of protons)
    //This solver follows the algo described in Thijssen, p74-76
    Bs = BS;
    nStates = Bs.Nstates;        //set number of states in basis
    nElectrons = N;              //set number of electrons
    nProtons = Z;                //set number of protons, may be removed

    //initializing all matrices and vectors
    //C.zeros(nStates,nElectrons/2);//set initial C equal to the unit matrix
    C.zeros(nStates,nStates); // Alternating between these two, alteration made to implement CoupledCluster extension of code
    Cu.zeros(nStates,nStates); // Alternating between these two, alteration made to implement CoupledCluster extension of code
    Cd.zeros(nStates,nStates); // Alternating between these two, alteration made to implement CoupledCluster extension of code

    //Cu.zeros(nStates,nElectronsU); // Alternating between these two, alteration made to implement CoupledCluster extension of code
    //Cd.zeros(nStates,nElectronsD); // Alternating between these two, alteration made to implement CoupledCluster extension of code


    F.zeros(nStates,nStates);     //initialize Fock matrix
    Fu.zeros(nStates,nStates);     //initialize Fock matrix
    Fd.zeros(nStates,nStates);     //initialize Fock matrix

    P.zeros(nStates,nStates);     //initialize Density matrix
    Pu.zeros(nStates,nStates);     //initialize Density matrix
    Pd.zeros(nStates,nStates);     //initialize Density matrix

    U.zeros(nStates,nStates);     //initialize Unitary matrix
    G.zeros(nStates,nStates);     //initialize Fock-component matrix

    Fprime.zeros(nStates,nStates);//transformed Fock matrix
    Fprimeu.zeros(nStates,nStates);//transformed Fock matrix
    Fprimed.zeros(nStates,nStates);//transformed Fock matrix

    epsilon.zeros(nStates);       //eigenvalues from diagonalization
    epsilon_prev.zeros(nStates);  //eigenvalues from previous diagonalization

    epsilonu.zeros(nStates);       //eigenvalues from diagonalization
    epsilonu_prev.zeros(nStates);  //eigenvalues from previous diagonalization

    epsilond.zeros(nStates);       //eigenvalues from diagonalization
    epsilond_prev.zeros(nStates);  //eigenvalues from previous diagonalization


    setupCoupledMatrix();  //import particle-particle interaction integrals
    setupP();              //set up initial density matrix
    s_diag.zeros(nStates);
}

double uhfsolve::solve(){
    //carefully following the steps laid out on pages 74-77 in Thijssen

    setupUnitMatrices();
    setupP();
    iterations = 0;

    setupF();
    //printMatrices();
    while(convergenceCriteria()){
        epsilon_prev = epsilon;
        energyPrev = energyCalc();

        setupF();

        diagonalizeF();
        normalizeC();
        updateP();

        iterations += 1;
        //cout << energy() << endl;
    }
    //cout << "Converged in " << iterations << " iterations." << endl;
    //printMatrices();
    //createDensityMap();
    //return iterations;
    //C = Cu + Cd;
    return energyCalc();
}

void uhfsolve::setupTotalCoefficientMatrix(){

}

double uhfsolve::evaluateProbabilityDensity(vec3 r){
    double result = 0;
    for(int p=0;p<nStates;p++){
        for(int q=0;q<nStates;q++){
            result += P(p,q)*Bs.evaluateContracted(p, r)*Bs.evaluateContracted(q,r);
        }
    }
    return result;
}

void uhfsolve::createDensityMap(string filename){
    int dim=1000;
    double dx = 0.01;
    double dV = dx*dx*dx;
    densityMap.zeros(dim,dim,dim);
    mat densitySlice;
    densitySlice.zeros(dim,dim);
    vec3 R;
    for(int i=0;i<dim;i++){
        for(int j=0;j<dim;j++){
            R = {(double) i*dx,(double) j*dx, 0.0};
            densitySlice(i,j)=evaluateProbabilityDensity(R);
            /*
            for(int k=-dim;k<dim;k++){
                R = {(double) i*dx,(double) j*dx, (double) k*dx};
                densityMap(i+dim,j+dim,k+dim) = evaluateProbabilityDensity(R)*dV;
            */
        }
    }
    //densityMap.print();
    //densityMap.save("testmap3", raw_ascii);
    //densitySlice.print();
    densitySlice.save(filename, raw_ascii);

}

void uhfsolve::setupCoupledMatrix(){
    //import particle-particle interaction integrals
    int n = Bs.Nstates;
    coupledMatrix.set_size(n, n);
    for (int p = 0; p<n; p++){
        for (int q = 0; q<n; q++){
            coupledMatrix(p, q) = zeros(n, n);
        }
    }
    for (int p = 0; p<n; p++){
        for (int q = 0; q<n; q++){
            for (int r = 0; r<n; r++){
                for (int s = 0; s<n; s++){
                    coupledMatrix(p, r)(q, s) = Bs.v(p, q)(r, s); //alt (1) "Strange" (CHANGED TODAY)

                }
            }
        }
    }
}

void uhfsolve::setupF(){
    //set up the Fock matrix
    double Di = 0;
    double Ex = 0;
    double DiMinusEx = 0;
    for(int p=0;p<nStates;p++){
        for(int q=0;q<nStates;q++){
            //F(p,q) = Bs.h(p,q);
            Fu(p,q) = Bs.h(p,q);
            Fd(p,q) = Bs.h(p,q);
            for(int r=0;r<nStates;r++){
                for(int s=0;s<nStates;s++){
                    //Di = Bs.v(p,r)(q,s);
                    //Ex = Bs.v(p,r)(s,q);

                    Di = Bs.v(p,q)(r,s);
                    Ex = Bs.v(p,s)(r,q);
                    DiMinusEx = Di-Ex;
                    //F(p,q) += 0.5*coupledMatrixTilde(p,q,r,s)*P(r,s);  //Alt. 2 27/5 2014
                    Fu(p,q) += Pu(s,r)*(Di-Ex) + Pd(s,r)*Di;
                    Fd(p,q) += Pd(s,r)*(Di-Ex) + Pu(s,r)*Di;
                }
            }
        }
    }
}

double uhfsolve::energyCalc(){
    //return energy for current Fock matrix

    //return 0.5*accu(P % (Bs.h + F))+Bs.nnInteraction();

    return 0.5*accu((Pu+Pd) % Bs.h + Fu % Pu + Fd % Pd)+Bs.nnInteraction();

    //return 0.5*accu((Pu.submat(0,0,nStates, nElectronsU)+Pd.submat(0,0,nStates,nElectronsD)) % Bs.h + Fu.submat(0,0,nStates, nElectronsU) % Pu.submat(0,0,nStates, nElectronsU) + Fd.submat(0,0,nStates,nElectronsD) % Pd.submat(0,0,nStates,nElectronsD))+Bs.nnInteraction();
}

double uhfsolve::energy(){
    //return ground state energy
    double e0 = 0;
    for(int p = 0; p<nStates;p++){
        for(int q = 0; q<nStates; q++){
            e0 += P(p,q)*Bs.h(p,q);
            for(int r = 0; r<nStates;r++){
                for(int s = 0; s<nStates; s++){
                    e0 += 0.25*coupledMatrixTilde(p,q,r,s)*P(p,q)*P(s,r);
                }
            }
        }
    }
    return e0+Bs.nnInteraction();
}

double uhfsolve::coupledMatrixTilde(int p, int q, int r, int s){
    //return direct and exchange term, weigthed to include spin
    //return 2*coupledMatrix(p,r)(q,s) - coupledMatrix(p,r)(s,q); //CHANGED 27/5
    return Bs.v(p,q)(r,s) - Bs.v(p,s)(r,q); //CHANGED 28/9
}

void uhfsolve::setupUnitMatrices(){
    //Bring overlap matrix to unit form, set up V
    eig_sym(s_diag,U,Bs.S);           //following Thijssen, p38-39
    V = U*diagmat(1.0/sqrt(s_diag));
}



void uhfsolve::diagonalizeF(){
    //diagonalize the Fock matrix
    Fprime = V.t()*F*V;
    Fprimeu = V.t()*Fu*V;
    Fprimed = V.t()*Fd*V;

    eig_sym(epsilon, Cprime, Fprime);
    eig_sym(epsilonu, Cprimeu, Fprimeu);
    eig_sym(epsilond, Cprimed, Fprimed);

    //C = V*Cprime.submat(0, 0, nStates - 1, nElectrons/2 -1);
    C = V*Cprime; //alternating between these two to implement coupled cluster calculations (need virtual orbitals)
    Cu = V*Cprimeu; //alternating between these two to implement coupled cluster calculations (need virtual orbitals)
    Cd = V*Cprimed; //alternating between these two to implement coupled cluster calculations (need virtual orbitals)
}

void uhfsolve::normalizeC(){
    //Normalize the coefficient matrix
    double norm, normu, normd;
    //for(int i = 0; i<nElectrons/2;i++){
    //    norm = dot(C.col(i),Bs.S*C.col(i));
    //    C.col(i) = C.col(i)/sqrt(norm);
    //}

    //Again, alternating the above with the following, as a preparation for the coupled cluster calculations:
    for(int i = 0; i<nStates;i++){
        norm = dot(C.col(i),Bs.S*C.col(i));
        C.col(i) = C.col(i)/sqrt(norm);

        normu = dot(Cu.col(i),Bs.S*Cu.col(i));
        Cu.col(i) = Cu.col(i)/sqrt(normu);

        normd = dot(Cd.col(i),Bs.S*Cd.col(i));
        Cd.col(i) = Cd.col(i)/sqrt(normd);
    }
}

void uhfsolve::setupP(){
    //setup density matrix, make a first guess
    //P.zeros(); //we don't have any reason to do this any other ways yet.
    //P.eye();
    //Pu.eye();
    //Pd.eye();
    //Pu.zeros(nStates, nElectronsU);
    //Pd.zeros(nStates, nElectronsD);


    //Pu = dampingFactor*Pu + (1-dampingFactor)*Cu.submat(0,0,nStates-1, nElectronsU-1)*Cu.submat(0,0,nStates-1, nElectronsU-1).t();
    //Pd = dampingFactor*Pd + (1-dampingFactor)*Cd.submat(0,0,nStates-1, nElectronsD-1)*Cd.submat(0,0,nStates-1, nElectronsD-1).t();

    //Pu = dampingFactor*Pu + (1-dampingFactor)*Cu.submat(0,0,nStates, nElectronsU)*Cu.submat(0,0,nStates, nElectronsU).t();
    //Pd = dampingFactor*Pd + (1-dampingFactor)*Cd.submat(0,0,nStates, nElectronsD)*Cd.submat(0,0,nStates, nElectronsD).t();

    //Pu = Cu.submat(0,0,nStates, nElectronsU)*Cu.submat(0,0,nStates, nElectronsU).t();
    //Pd = Cd.submat(0,0,nStates, nElectronsD)*Cd.submat(0,0,nStates, nElectronsD).t();;

    //Pu = Cu*Cu.t();
    //Pd = Cd*Cd.t();

    Pu = Cu.submat(0,0,nStates-1, nElectronsU-1)*Cu.submat(0,0,nStates-1, nElectronsU-1).t(); //new
    Pu(0,1) = 0.1; //non-interacting initial condition
    Pd = Cd.submat(0,0,nStates-1, nElectronsD-1)*Cd.submat(0,0,nStates-1, nElectronsD-1).t();
}

void uhfsolve::updateP(){
    //construct new density matrix
    P = dampingFactor*P + (1-dampingFactor)*2.0*C.cols(0, nElectrons/2.0 -1)*C.cols(0, nElectrons/2.0 -1).t();

    //Pu = dampingFactor*Pu + (1-dampingFactor)*Cu*Cu.t(); //will produce correct results for nElectronsU!=nElectronsD
    //Pd = dampingFactor*Pd + (1-dampingFactor)*Cd*Cd.t();

    //Pu = dampingFactor*Pu + (1-dampingFactor)*Cu.submat(0,0,nStates, nElectronsU)*Cu.submat(0,0,nStates, nElectronsU).t(); //will produce correct result for spin-symmetric systems
    //Pd = dampingFactor*Pd + (1-dampingFactor)*Cd.submat(0,0,nStates, nElectronsD)*Cd.submat(0,0,nStates, nElectronsD).t();

    Pu = dampingFactor*Pu + (1-dampingFactor)*Cu.submat(0,0,nStates-1, nElectronsU-1)*Cu.submat(0,0,nStates-1, nElectronsU-1).t(); //new
    Pd = dampingFactor*Pd + (1-dampingFactor)*Cd.submat(0,0,nStates-1, nElectronsD-1)*Cd.submat(0,0,nStates-1, nElectronsD-1).t();
}

bool uhfsolve::convergenceCriteria(){
    //Evaluate convergence conditions
    bool condition = true;
    if(iterations>5000){
        condition = false;
    }
    if(abs(energyPrev-energy())<tolerance){

        condition = false;
    }
    return condition;
}


void uhfsolve::setupCoupledMatrix_unused(){
    //still following Dragly, further references to indexation differences between Thijssen and Helgaker
    coupledMatrix.set_size(nStates,nStates);
    for(int p = 0; p<nStates;p++){
        for(int r = 0; r<nStates;r++){
            coupledMatrix(p,r)=zeros(nStates,nStates);
            for(int q = p; q<nStates;q++){
                for(int s = r; s<nStates;s++){
                    coupledMatrix(p,r)(q,s) = Bs.v(p,q)(r,s);
                }
            }
        }
    }
    double val;
    for(int p = 0; p<nStates;p++){
        for(int r = 0; r<nStates;r++){
            for(int q = p; q<nStates;q++){
                for(int s = r; s<nStates;s++){
                    val = coupledMatrix(p,r)(q,s);
                    coupledMatrix(q,s)(p,r) = val;
                    coupledMatrix(q,r)(p,s) = val;
                    coupledMatrix(p,s)(q,r) = val;
                    coupledMatrix(r,p)(s,q) = val;
                    coupledMatrix(s,p)(r,q) = val;
                    coupledMatrix(r,q)(s,p) = val;
                    coupledMatrix(s,q)(r,p) = val;
                }
            }
        }
    }
    coupledMatrix.print();
}

void uhfsolve::printMatrices(){
    cout << "Fock matrix" << endl;
    F.print();
    cout << " " << endl;

    cout << "Coeff matrix" << endl;
    C.print();
    cout << " " << endl;

    cout << "Density matrix" << endl;
    P.print();
    cout << " " << endl;
    cout << "----------------------" << endl;

    cout << "H matrix" << endl;
    Bs.h.print();
    cout << " " << endl;
    cout << "----------------------" << endl;


    cout << "Coupled matrix" << endl;
    coupledMatrix.print();
    cout << " " << endl;
    cout << "----------------------" << endl;

    cout << "Overlap" << endl;
    Bs.S.print();
    cout << " " << endl;
    cout << "----------------------" << endl;



}

