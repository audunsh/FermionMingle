#include "basis/basis.h"
#include "basis/contracted.h"
#include "basis/primitive.h"
#include "basis/electrongas.h"
#include <armadillo>
#include "integrator/integrator.h"
#include "integrator/boysfunction.h"


basis::basis(){
    set_size(3);
}


void basis::setup_electrongas(int NEv, double Lv){
    L = Lv;
    NE = NEv;
    gasbasis.generate_state_list(NE, L);
    set_size(gasbasis.n_basis_functions);
    //Set up the interactions for an electron gas with NE electrons in a cube with volume L^3
}

void basis::init_electron_integrals(){
    //Calculate all integrals for electron interactions (needed only for testing)
    for(int p=0; p<Nstates; p++){
        for(int q=0; q<Nstates; q++){

            S(p,q) = 1*(p==q);
            h(p,q) = gasbasis.h(p,q);
            for(int r=0; r<Nstates; r++){
                for(int s=0; s<Nstates; s++){
                    v(p,q)(r,s) = gasbasis.v(p,q,r,s);
                }
            }
        }
        //cout << p << endl;
    }
    //h.print();
    //v.print();
}


void basis::set_orthonormal(){
    //sets the overlap matrix equal to the identity matrix
    for(int i=0;i < Nstates;i++){
        for(int j=0;j < Nstates;j++){
            if(i==j){
                S(i,j) = 1.0;
            }
        }
    }
}


void basis::read(string filename, int Zn){
    Z = Zn;
    //read basis from file
    //this also means you have to provide the one-body energy and the overlap matrix (use h0 and ...)
    cout << "Loading predefined basis from file: " << filename << endl;
    ifstream myfile;
    myfile.open(filename.c_str());
    if(!myfile.good()){
        cout << "Failed to load file." << endl;
    }
    if (myfile.is_open()){
        int p,q,r,s;
        double value;
        while (!myfile.eof()){
            myfile >> p;
            myfile >> q;
            myfile >> r;
            myfile >> s;
            myfile >> value;
            try{
                v(p,q)(r,s) = Z*value;
            }
            catch(int e){
                cout << "Failed to read basis from file." << endl;
            }
        }
    }
    else
        cout << "Did not manage to load basis from file."<< endl;
    for(int p=0; p<Nstates;p++){
        for(int q=0; q<Nstates;q++){
            h(p,q) = h0(p,q);
        }
    }
    set_orthonormal();
}

void basis::expand(){
    //expand basis for explicid spin inclusion
    Nstates2 = 2*Nstates;
    field<mat> Vn(Nstates2, Nstates2); //spin-basis
    for (int i = 0; i < Nstates2; ++i) {
        for (int j = 0; j < Nstates2; ++j) {
            Vn(i,j) = zeros<mat>(Nstates2,Nstates2); // fill V with 3x3 mx elements
        }
    }
    V = Vn;
    h.set_size(Nstates2);
    S.set_size(Nstates2,Nstates2);
    h.zeros();
    S.zeros();
    H.set_size(Nstates2,Nstates2);
    H.zeros();
    double D = 0;
    double Ex = 0;
    for (int p = 0; p < Nstates2; ++p) {
        for (int q = 0; q < Nstates2; ++q) {
            for (int r = 0; r < Nstates2; ++r) {
                for (int s= 0; s < Nstates2; ++s) {
                    try{
                        D = v(p/2,q/2)(r/2,s/2);  // Direct term
                        Ex =v(p/2,q/2)(s/2,r/2);  // Exchange term
                        V(p,q)(r,s) = state(p,q,r,s,D,Ex);
                    }
                    catch(int e){
                        cout << "Failed to read basis." << endl;
                    }
                }
            }
            H(p,q) = h0(p,q);
        }
    }

}

double basis::h0(int i, int j){
    // the one-body interaction
    double h = 0;
    if (i == j){
        double n = i + 1.0;
        h = -(Z*Z)/(2*n*n);
    }
    return h;
}

double basis::state(int p, int q, int r, int s, double D, double Ex){
    //Explicit spin implementation using direct and exchange term
    double S = 0;
    int s1 = 0; int s2 = 0;int s3 = 0;int s4 = 0;
    s1 = p%2;
    s2 = q%2;
    s3 = r%2;
    s4 = s%2;
    if (s1 == s2){
        if (s3 == s4){
            if ( s1 == s3){
                S = D-Ex;
            }
            else{
                S = 0;
            }
        }
    }
    else if (s1 != s2){
        if (s3 != s4){
            if (s1 == s3){
                S = D;
            }
            else{
                S = -Ex;
            }
        }
        else{
            S = 0;
        }
    }
    return S;
}
void basis::add_state(){
    Nstates += 1;
    nPrimitivesInState.set_size(Nstates);
    basisSts.push_back(contracted());
}

void basis::add_primitive_to_state(int Stateindex, Primitive P){
    basisSts.at(Stateindex).appendPrimitive(P);
    nPrimitivesInState(Stateindex) += 1;
}

void basis::add_nucleus(vec3 pos, int charge){
    int Nsize = nucleusCharges.size();
    nucleusCharges.set_size(Nsize + 1);
    nucleusPositions.set_size(Nsize + 1);
    nucleusPositions(Nsize) = pos;
    nucleusCharges(Nsize) = charge;
}

void basis::add_STO_NG_s_orbital(int NG, vec weigths, vec exponents, vec3 corePos){
    add_state();
    for(int i = 0; i<NG;i++){
        Primitive S = turbomolePrimitive(weigths(i),exponents(i),0,0,0,corePos);
        add_primitive_to_state(Nstates-1, S);
    }
}

void basis::add_STO_NG_p_orbital(int NG, vec exponents, vec weigths, vec3 corePos){
    add_state();
    for(int i = 0; i<NG;i++){
        Primitive P1 = turbomolePrimitive(weigths(i),exponents(i),1,0,0,corePos);
        add_primitive_to_state(Nstates-1, P1);
    }
    add_state();
    for(int i = 0; i<NG;i++){
        Primitive P2 = turbomolePrimitive(weigths(i),exponents(i),0,1,0,corePos);
        add_primitive_to_state(Nstates-1, P2);
    }
    add_state();
    for(int i = 0; i<NG;i++){
        Primitive P3 = turbomolePrimitive(weigths(i),exponents(i),0,0,1,corePos);
        add_primitive_to_state(Nstates-1, P3);
    }
}

void basis::add_STO_NG_d_orbital(int NG, vec exponents, vec weigths, vec3 corePos){
    add_state();
    for(int i = 0; i<NG;i++){
        Primitive D1 = turbomolePrimitive(weigths(i),exponents(i),2,0,0,corePos);
        add_primitive_to_state(Nstates-1, D1);
    }

    add_state();
    for(int i = 0; i<NG;i++){
        Primitive D2 = turbomolePrimitive(weigths(i),exponents(i),0,2,0,corePos);
        add_primitive_to_state(Nstates-1, D2);
    }

    add_state();
    for(int i = 0; i<NG;i++){
        Primitive D3 = turbomolePrimitive(weigths(i),exponents(i),0,0,2,corePos);
        add_primitive_to_state(Nstates-1, D3);
    }

    add_state();
    for(int i = 0; i<NG;i++){
        Primitive D4 = turbomolePrimitive(weigths(i),exponents(i),1,1,0,corePos);
        add_primitive_to_state(Nstates-1, D4);
    }

    add_state();
    for(int i = 0; i<NG;i++){
        Primitive D5 = turbomolePrimitive(weigths(i),exponents(i),0,1,1,corePos);
        add_primitive_to_state(Nstates-1, D5);
    }

    add_state();
    for(int i = 0; i<NG;i++){
        Primitive D6 = turbomolePrimitive(weigths(i),exponents(i),1,0,1,corePos);
        add_primitive_to_state(Nstates-1, D6);
    }
}


void basis::init_Be2(vec3 corePos1, vec3 corePos2){
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(2);
    nucleusCharges.set_size(2);
    nucleusPositions(0) = corePos1;
    nucleusPositions(1) = corePos2;
    nucleusCharges(0) = 4;
    nucleusCharges(1) = 4;
    add_atom_STO3G("Be", corePos1);
    add_atom_STO3G("Be", corePos2);
    set_size(Nstates);
}

void basis::init_H2(vec3 corePos1, vec3 corePos2){
    Nstates = 0;
    Nprimitives = 3;

    //nucleusPositions.set_size(0);
    //nucleusCharges.set_size(0);
    //nucleusPositions(0) = corePos1;
    //nucleusPositions(1) = corePos2;
    //nucleusCharges(0) = 1;
    //nucleusCharges(1) = 1;

    add_nucleus(corePos1, 1);
    add_nucleus(corePos2, 1);

    add_atom_STO3G("H", corePos1);
    add_atom_STO3G("H", corePos2);

    set_size(Nstates);

}

void basis::init_He2(vec3 corePos1, vec3 corePos2){
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(2);
    nucleusCharges.set_size(2);
    nucleusPositions(0) = corePos1;
    nucleusPositions(1) = corePos2;
    nucleusCharges(0) = 2;
    nucleusCharges(1) = 2;
    add_atom_STO3G("He", corePos1);
    add_atom_STO3G("He", corePos2);
    set_size(Nstates);
}

void basis::init_H2O(vec3 corePos1, vec3 corePos2, vec3 corePos3){
    //basisSts.clear();
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(3);
    nucleusCharges.set_size(3);
    nucleusPositions(0) = corePos1;
    nucleusPositions(1) = corePos2;
    nucleusPositions(2) = corePos3;
    nucleusCharges(0) = 1;
    nucleusCharges(1) = 1;
    nucleusCharges(2) = 8;
    add_atom_STO3G("H", corePos1);
    add_atom_STO3G("H", corePos2);
    add_atom_STO3G("O", corePos3);
    set_size(Nstates);
}

void basis::init_H2Be(vec3 corePos1, vec3 corePos2, vec3 corePos3){
    //basisSts.clear();
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(3);
    nucleusCharges.set_size(3);
    nucleusPositions(0) = corePos1;
    nucleusPositions(1) = corePos2;
    nucleusPositions(2) = corePos3;
    nucleusCharges(0) = 1;
    nucleusCharges(1) = 1;
    nucleusCharges(2) = 4;
    add_atom_STO3G("H", corePos1);
    add_atom_STO3G("H", corePos2);
    add_atom_STO3G("Be", corePos3);
    set_size(Nstates);
}


void basis::init_O2(vec3 corePos1, vec3 corePos2){
    //basisSts.clear();
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(2);
    nucleusCharges.set_size(2);
    nucleusPositions(0) = corePos1;
    nucleusPositions(1) = corePos2;
    nucleusCharges(0) = 8;
    nucleusCharges(1) = 8;
    add_atom_STO3G("O", corePos1);
    add_atom_STO3G("O", corePos2);
    set_size(Nstates);
}

void basis::init_O(vec3 corePos1){
    //basisSts.clear();
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(1);
    nucleusCharges.set_size(1);
    nucleusPositions(0) = corePos1;
    nucleusCharges(0) = 8.0;
    add_atom_STO3G("O", corePos1);
    set_size(Nstates);
}

void basis::init_Ne(vec3 corePos1){
    //basisSts.clear();
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(1);
    nucleusCharges.set_size(1);
    nucleusPositions(0) = corePos1;
    nucleusCharges(0) = 10;
    add_atom_STO3G("Ne", corePos1);
    set_size(Nstates);
}

void basis::init_He(vec3 corePos1){
    //basisSts.clear();
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(1);
    nucleusCharges.set_size(1);
    nucleusPositions(0) = corePos1;
    nucleusCharges(0) = 2.0;
    add_atom_STO3G("He", corePos1);
    set_size(Nstates);
}

void basis::init_Be(vec3 corePos1){
    //basisSts.clear();
    Nstates = 0;
    Nprimitives = 3;
    nucleusPositions.set_size(1);
    nucleusCharges.set_size(1);
    nucleusPositions(0) = corePos1;
    nucleusCharges(0) = 4.0;
    add_atom_STO3G("Be", corePos1);
    set_size(Nstates);
}


Primitive basis::turbomolePrimitive(double weight, double exponent, double i, double j, double k, vec3 corePos){
    //return Primitive(turboNormalization(exponent,i,j,k)*weight,i,j,k,exponent,corePos);
    return Primitive(turboNormalization(exponent,i,j,k)*weight,i,j,k,exponent,corePos);

}

//void basis::import(TurboMoleParser system, vec3 corePos){
    //imports one state basis from a file using a parser written by Svenn-Arne Dragly
//}

void basis::add_atom_STO3G(string configuration, vec3 corePos){
    //add new atom to the basis
    //turbo(we,exp)
    if(configuration == "H"){
        vec weigths =   {0.15432897, 0.53532814,0.44463454};
        vec exponents   = {3.42525091, 0.62391373,0.16885540};
        add_STO_NG_s_orbital(3 , weigths,exponents,corePos);
    }

    if(configuration == "He"){
        vec weigths =   {0.15432897,0.53532814,0.44463454};
        vec exponents   = {6.36242139,1.15892300,0.31364979};
        add_STO_NG_s_orbital(3 , weigths,exponents,corePos);
    }

    if(configuration == "Be"){
        add_state();
        Primitive S1A = turbomolePrimitive(0.15432897,30.1678710,0,0,0,corePos);
        Primitive S1B = turbomolePrimitive(0.53532814,5.4951153,0,0,0,corePos);
        Primitive S1C = turbomolePrimitive(0.44463454,1.4871927,0,0,0,corePos);
        add_primitive_to_state(Nstates-1, S1A);
        add_primitive_to_state(Nstates-1, S1B);
        add_primitive_to_state(Nstates-1, S1C);

        add_state();
        Primitive S2A = turbomolePrimitive(-0.09996723,1.3148331,0,0,0,corePos);
        Primitive S2B = turbomolePrimitive(0.39951283,0.3055389,0,0,0,corePos);
        Primitive S2C = turbomolePrimitive(0.70011547,0.0993707,0,0,0,corePos);
        add_primitive_to_state(Nstates-1, S2A);
        add_primitive_to_state(Nstates-1, S2B);
        add_primitive_to_state(Nstates-1, S2C);

        add_state();
        Primitive P1A = turbomolePrimitive(0.15591627,1.3148331,1,0,0,corePos);
        Primitive P1B = turbomolePrimitive(0.60768372,0.3055389,1,0,0,corePos);
        Primitive P1C = turbomolePrimitive(0.39195739,0.0993707,1,0,0,corePos);
        add_primitive_to_state(Nstates-1, P1A);
        add_primitive_to_state(Nstates-1, P1B);
        add_primitive_to_state(Nstates-1, P1C);

        add_state();
        Primitive P2A = turbomolePrimitive(0.15591627,1.3148331,0,1,0,corePos);
        Primitive P2B = turbomolePrimitive(0.60768372,0.3055389,0,1,0,corePos);
        Primitive P2C = turbomolePrimitive(0.39195739,0.0993707,0,1,0,corePos);
        add_primitive_to_state(Nstates-1, P2A);
        add_primitive_to_state(Nstates-1, P2B);
        add_primitive_to_state(Nstates-1, P2C);

        add_state();
        Primitive P3A = turbomolePrimitive(0.15591627,1.3148331,0,0,1,corePos);
        Primitive P3B = turbomolePrimitive(0.60768372,0.3055389,0,0,1,corePos);
        Primitive P3C = turbomolePrimitive(0.39195739,0.0993707,0,0,1,corePos);
        add_primitive_to_state(Nstates-1, P3A);
        add_primitive_to_state(Nstates-1, P3B);
        add_primitive_to_state(Nstates-1, P3C);

    }

    if(configuration == "Ne"){
        add_state();
        Primitive S1A = turbomolePrimitive(0.15432897,207.0156100,0,0,0,corePos);
        Primitive S1B = turbomolePrimitive(0.53532814,37.7081510,0,0,0,corePos);
        Primitive S1C = turbomolePrimitive(0.44463454,10.2052970,0,0,0,corePos);
        add_primitive_to_state(Nstates-1, S1A);
        add_primitive_to_state(Nstates-1, S1B);
        add_primitive_to_state(Nstates-1, S1C);


        add_state();
        Primitive S2A = turbomolePrimitive(-0.09996723,8.2463151,0,0,0,corePos);
        Primitive S2B = turbomolePrimitive(0.39951283,1.9162662,0,0,0,corePos);
        Primitive S2C = turbomolePrimitive(0.70011547,0.6232293,0,0,0,corePos);
        add_primitive_to_state(Nstates-1, S2A);
        add_primitive_to_state(Nstates-1, S2B);
        add_primitive_to_state(Nstates-1, S2C);

        add_state();
        Primitive P1A = turbomolePrimitive(0.15591627,8.2463151,1,0,0,corePos);
        Primitive P1B = turbomolePrimitive(0.60768372,1.9162662,1,0,0,corePos);
        Primitive P1C = turbomolePrimitive(0.39195739,0.6232293,1,0,0,corePos);
        add_primitive_to_state(Nstates-1, P1A);
        add_primitive_to_state(Nstates-1, P1B);
        add_primitive_to_state(Nstates-1, P1C);

        add_state();
        Primitive P2A = turbomolePrimitive(0.15591627,8.2463151,0,1,0,corePos);
        Primitive P2B = turbomolePrimitive(0.60768372,1.9162662,0,1,0,corePos);
        Primitive P2C = turbomolePrimitive(0.39195739,0.6232293,0,1,0,corePos);
        add_primitive_to_state(Nstates-1, P2A);
        add_primitive_to_state(Nstates-1, P2B);
        add_primitive_to_state(Nstates-1, P2C);

        add_state();
        Primitive P3A = turbomolePrimitive(0.15591627,8.2463151,0,0,1,corePos);
        Primitive P3B = turbomolePrimitive(0.60768372,1.9162662,0,0,1,corePos);
        Primitive P3C = turbomolePrimitive(0.39195739,0.6232293,0,0,1,corePos);
        add_primitive_to_state(Nstates-1, P3A);
        add_primitive_to_state(Nstates-1, P3B);
        add_primitive_to_state(Nstates-1, P3C);
    }

    if(configuration == "O"){
        add_state();
        Primitive S1A = turbomolePrimitive(0.15432897,130.7093200,0,0,0,corePos);
        Primitive S1B = turbomolePrimitive(0.53532814,23.8088610,0,0,0,corePos);
        Primitive S1C = turbomolePrimitive(0.44463454,6.4436083,0,0,0,corePos);
        add_primitive_to_state(Nstates-1, S1A);
        add_primitive_to_state(Nstates-1, S1B);
        add_primitive_to_state(Nstates-1, S1C);

        add_state();
        Primitive S2A = turbomolePrimitive(-0.09996723,5.0331513,0,0,0,corePos);
        Primitive S2B = turbomolePrimitive(0.39951283,1.1695961,0,0,0,corePos);
        Primitive S2C = turbomolePrimitive(0.70011547,0.3803890,0,0,0,corePos);
        add_primitive_to_state(Nstates-1, S2A);
        add_primitive_to_state(Nstates-1, S2B);
        add_primitive_to_state(Nstates-1, S2C);

        add_state();
        Primitive P1A = turbomolePrimitive(0.15591627,5.0331513,1,0,0,corePos);
        Primitive P1B = turbomolePrimitive(0.60768372,1.1695961,1,0,0,corePos);
        Primitive P1C = turbomolePrimitive(0.39195739,0.3803890,1,0,0,corePos);
        add_primitive_to_state(Nstates-1, P1A);
        add_primitive_to_state(Nstates-1, P1B);
        add_primitive_to_state(Nstates-1, P1C);

        add_state();
        Primitive P2A = turbomolePrimitive(0.15591627,5.0331513,0,1,0,corePos);
        Primitive P2B = turbomolePrimitive(0.60768372,1.1695961,0,1,0,corePos);
        Primitive P2C = turbomolePrimitive(0.39195739,0.3803890,0,1,0,corePos);
        add_primitive_to_state(Nstates-1, P2A);
        add_primitive_to_state(Nstates-1, P2B);
        add_primitive_to_state(Nstates-1, P2C);

        add_state();
        Primitive P3A = turbomolePrimitive(0.15591627,5.0331513,0,0,1,corePos);
        Primitive P3B = turbomolePrimitive(0.60768372,1.1695961,0,0,1,corePos);
        Primitive P3C = turbomolePrimitive(0.39195739,0.3803890,0,0,1,corePos);
        add_primitive_to_state(Nstates-1, P3A);
        add_primitive_to_state(Nstates-1, P3B);
        add_primitive_to_state(Nstates-1, P3C);
    }

    if(configuration == "Si"){
        Nstates += 5;
        Primitive S1A(turboNormalization(407.7975514,0,0,0)*0.1543289673,0,0,0,407.7975514,corePos);
        Primitive S1B(turboNormalization(74.28083305,0,0,0)*0.5353281423,0,0,0,74.28083305, corePos);
        Primitive S1C(turboNormalization(20.10329229,0,0,0)*0.4446345422,0,0,0,20.10329229, corePos);

        Primitive S2A(turboNormalization(23.19365606,0,0,0)*-0.09996722919,0,0,0,23.19365606,corePos);
        Primitive S2B(turboNormalization(5.389706871,0,0,0)*0.39951282610,0,0,0,5.389706871, corePos);
        Primitive S2C(turboNormalization(1.752899952,0,0,0)*0.70011546890,0,0,0,1.752899952, corePos);

        Primitive S3A(turboNormalization(1.4787406220,0,0,0)*-0.2196203690,0,0,0,1.4787406220,corePos);
        Primitive S3B(turboNormalization(0.4125648801,0,0,0)*0.2255954336,0,0,0,0.4125648801, corePos);
        Primitive S3C(turboNormalization(1.752899952,0,0,0)*0.9003984260,0,0,0,1.752899952, corePos);

        Primitive P1A(turboNormalization(23.19365606,1,0,0)*0.1559162750,1,0,0,23.19365606, corePos);
        Primitive P1B(turboNormalization(5.389706871,1,0,0)*0.6076837186,1,0,0,5.389706871, corePos);
        Primitive P1C(turboNormalization(1.752899952,1,0,0)*0.3919573931,1,0,0,1.752899952, corePos);



        Primitive P2A(turboNormalization(1.4787406220,0,1,0)*0.01058760429,0,1,0,1.4787406220, corePos);
        Primitive P2B(turboNormalization(0.4125648801,0,1,0)*0.59516700530,0,1,0,0.4125648801, corePos);
        Primitive P2C(turboNormalization(0.1614750979,0,1,0)*0.46200101200,0,1,0,0.1614750979, corePos);

        Primitive S1[3] = {S1A,S1B,S1C};
        Primitive S2[3] = {S2A,S2B,S2C};
        Primitive S3[3] = {S3A,S3B,S3C};

        Primitive P1[3] = {P1A,P1B,P1C};
        Primitive P2[3] = {P2A,P2B,P2C};

        contracted C1(3,S1);
        contracted C2(3,S2);
        contracted C3(3,S3);
        contracted C4(3,P1);
        contracted C5(3,P2);

        basisSts.push_back(C1);
        basisSts.push_back(C2);
        basisSts.push_back(C3);
        basisSts.push_back(C4);
        basisSts.push_back(C5);
    }
}

void basis::init_molecule(string configuration, vec nProtons, field<vec> corePos){
    //using STO-3G for molecules
    if(configuration == "Be"){
        //multiple Be atoms interacting as a molecule
        Nstates = 3*nProtons.size();
        set_size(Nstates);
        Nprimitives = 3;
        nucleusPositions.set_size(nProtons.size());
        for(int i=0;i<nProtons.size();i++){
            nucleusCharges(i) = nProtons(i);
            nucleusPositions(i) = corePos(i);
            Primitive S1A(turboNormalization(30.1678710,0,0,0)*0.15432897,0,0,0,30.1678710,corePos(i));
            Primitive S1B(turboNormalization(5.4951153,0,0,0)*0.53532814,0,0,0,5.4951153, corePos(i));
            Primitive S1C(turboNormalization(1.4871927,0,0,0)*0.44463454,0,0,0,1.4871927, corePos(i));

            Primitive S2A(turboNormalization(1.3148331,0,0,0)*-0.09996723,0,0,0,1.3148331,corePos(i));
            Primitive S2B(turboNormalization(0.3055389,0,0,0)*0.39951283,0,0,0,0.3055389, corePos(i));
            Primitive S2C(turboNormalization(0.0993707,0,0,0)*0.70011547,0,0,0,0.0993707, corePos(i));

            Primitive P1A(turboNormalization(1.3148331,1,0,0)*0.15591627,1,0,0,1.3148331, corePos(i));
            Primitive P1B(turboNormalization(0.3055389,1,0,0)*0.60768372,1,0,0,0.3055389, corePos(i));
            Primitive P1C(turboNormalization(0.0993707,1,0,0)*0.39195739,1,0,0,0.0993707, corePos(i));

            Primitive S1[3] = {S1A,S1B,S1C};
            Primitive S2[3] = {S2A,S2B,S2C};
            Primitive P1[3] = {P1A,P1B,P1C};

            contracted C1 (3,S1);
            contracted C2(3,S2);
            contracted C3(3,P1);

            basisSts.push_back(C1);
            basisSts.push_back(C2);
            basisSts.push_back(C3);
        }
    }
    if(configuration == "H"){
        //multiple H-atoms interacting as a molecule
        Nstates = nProtons.size();
        set_size(Nstates);
        Nprimitives = 3;
        nucleusPositions.set_size(nProtons.size());
        for(int i=0;i<nProtons.size();i++){
            nucleusCharges(i) = nProtons(i);
            nucleusPositions(i) = corePos(i);
            Primitive S1A(turboNormalization(6.36242139,0,0,0)*0.15432897,0,0,0,6.36242139,corePos(i));
            Primitive S1B(turboNormalization(1.15892300,0,0,0)*0.53532814,0,0,0,1.15892300,corePos(i));
            Primitive S1C(turboNormalization(0.31364979,0,0,0)*0.44463454,0,0,0,0.31364979,corePos(i));

            Primitive S1[3] = {S1A,S1B,S1C};

            contracted C1 (3,S1);

            basisSts.push_back(C1);
        }
    }
    if(configuration == "O"){
        //multiple H-atoms interacting as a molecule
        Nstates = 3*nProtons.size();
        set_size(Nstates);
        Nprimitives = 3;
        nucleusPositions.set_size(nProtons.size());
        for(int i=0;i<nProtons.size();i++){
            nucleusCharges(i) = nProtons(i);
            nucleusPositions(i) = corePos(i);
            Primitive S1A(turboNormalization(130.7093200,0,0,0)*0.15432897,0,0,0,130.7093200,corePos(i));
            Primitive S1B(turboNormalization(23.8088610,0,0,0)*0.53532814,0,0,0,23.8088610, corePos(i));
            Primitive S1C(turboNormalization(6.4436083,0,0,0)*0.44463454,0,0,0,6.4436083, corePos(i));

            Primitive S2A(turboNormalization(5.0331513,0,0,0)*-0.09996723,0,0,0,5.0331513,corePos(i));
            Primitive S2B(turboNormalization(1.1695961,0,0,0)*0.39951283,0,0,0,1.1695961, corePos(i));
            Primitive S2C(turboNormalization(0.3803890,0,0,0)*0.70011547,0,0,0,0.3803890, corePos(i));

            Primitive P1A(turboNormalization(5.0331513,1,0,0)*0.15591627,1,0,0,5.0331513, corePos(i));
            Primitive P1B(turboNormalization(1.1695961,1,0,0)*0.60768372,1,0,0,1.1695961, corePos(i));
            Primitive P1C(turboNormalization(0.3803890,1,0,0)*0.39195739,1,0,0,0.3803890, corePos(i));

            Primitive S1[3] = {S1A,S1B,S1C};
            Primitive S2[3] = {S2A,S2B,S2C};
            Primitive P1[3] = {P1A,P1B,P1C};

            contracted C1 (3,S1);
            contracted C2(3,S2);
            contracted C3(3,P1);

            basisSts.push_back(C1);
            basisSts.push_back(C2);
            basisSts.push_back(C3);
        }
    }

    if(configuration == "Si"){
        //multiple Si-atoms interacting as a molecule
        Nstates = 5*nProtons.size();
        set_size(Nstates);
        Nprimitives = 3;
        nucleusPositions.set_size(nProtons.size());
        for(int i=0;i<nProtons.size();i++){
            nucleusCharges(i) = nProtons(i);
            nucleusPositions(i) = corePos(i);
            Primitive S1A(turboNormalization(407.7975514,0,0,0)*0.1543289673,0,0,0,407.7975514,corePos(i));
            Primitive S1B(turboNormalization(74.28083305,0,0,0)*0.5353281423,0,0,0,74.28083305, corePos(i));
            Primitive S1C(turboNormalization(20.10329229,0,0,0)*0.4446345422,0,0,0,20.10329229, corePos(i));

            Primitive S2A(turboNormalization(23.19365606,0,0,0)*-0.09996722919,0,0,0,23.19365606,corePos(i));
            Primitive S2B(turboNormalization(5.389706871,0,0,0)*0.39951282610,0,0,0,5.389706871, corePos(i));
            Primitive S2C(turboNormalization(1.752899952,0,0,0)*0.70011546890,0,0,0,1.752899952, corePos(i));

            Primitive S3A(turboNormalization(1.4787406220,0,0,0)*-0.2196203690,0,0,0,1.4787406220,corePos(i));
            Primitive S3B(turboNormalization(0.4125648801,0,0,0)*0.2255954336,0,0,0,0.4125648801, corePos(i));
            Primitive S3C(turboNormalization(1.752899952,0,0,0)*0.9003984260,0,0,0,1.752899952, corePos(i));

            Primitive P1A(turboNormalization(23.19365606,1,0,0)*0.1559162750,1,0,0,23.19365606, corePos(i));
            Primitive P1B(turboNormalization(5.389706871,1,0,0)*0.6076837186,1,0,0,5.389706871, corePos(i));
            Primitive P1C(turboNormalization(1.752899952,1,0,0)*0.3919573931,1,0,0,1.752899952, corePos(i));

            Primitive P2A(turboNormalization(1.4787406220,0,1,0)*0.01058760429,0,1,0,1.4787406220, corePos(i));
            Primitive P2B(turboNormalization(0.4125648801,0,1,0)*0.59516700530,0,1,0,0.4125648801, corePos(i));
            Primitive P2C(turboNormalization(0.1614750979,0,1,0)*0.46200101200,0,1,0,0.1614750979, corePos(i));

            Primitive S1[3] = {S1A,S1B,S1C};
            Primitive S2[3] = {S2A,S2B,S2C};
            Primitive S3[3] = {S3A,S3B,S3C};

            Primitive P1[3] = {P1A,P1B,P1C};
            Primitive P2[3] = {P2A,P2B,P2C};

            contracted C1(3,S1);
            contracted C2(3,S2);
            contracted C3(3,S3);
            contracted C4(3,P1);
            contracted C5(3,P2);

            basisSts.push_back(C1);
            basisSts.push_back(C2);
            basisSts.push_back(C3);
            basisSts.push_back(C4);
            basisSts.push_back(C5);

        }
    }
}

void basis::init_STO_3G(string configuration, double nProtons){
    //initialize STO-3G basis sets, following slides from Helgaker (2006)
    Z = nProtons; //set nuclear charge
    nucleusCharges.set_size(1);
    nucleusCharges(0) = nProtons;
    nucleusPositions.set_size(1);
    nucleusPositions(0) = {0,0,0};

    if(configuration == "Be"){
        //basisSet[3];
        Nstates = 5;
        set_size(Nstates);
        Nprimitives = 3;

        Primitive S1A(turboNormalization(30.1678710,0,0,0)*0.15432897,0,0,0,30.1678710,{0,0,0});
        Primitive S1B(turboNormalization(5.4951153,0,0,0)*0.53532814,0,0,0,5.4951153, {0,0,0});
        Primitive S1C(turboNormalization(1.4871927,0,0,0)*0.44463454,0,0,0,1.4871927, {0,0,0});


        Primitive S2A(turboNormalization(1.3148331,0,0,0)*-0.09996723,0,0,0,1.3148331,{0,0,0});
        Primitive S2B(turboNormalization(0.3055389,0,0,0)*0.39951283,0,0,0,0.3055389, {0,0,0});
        Primitive S2C(turboNormalization(0.0993707,0,0,0)*0.70011547,0,0,0,0.0993707, {0,0,0});

        Primitive P1A(turboNormalization(1.3148331,1,0,0)*0.15591627,1,0,0,1.3148331,{0,0,0});
        Primitive P1B(turboNormalization(0.3055389,1,0,0)*0.60768372,1,0,0,0.3055389,{0,0,0});
        Primitive P1C(turboNormalization(0.0993707,1,0,0)*0.39195739,1,0,0,0.0993707,{0,0,0});

        Primitive P2A(turboNormalization(1.3148331,0,1,0)*0.15591627,0,1,0,1.3148331,{0,0,0});
        Primitive P2B(turboNormalization(0.3055389,0,1,0)*0.60768372,0,1,0,0.3055389,{0,0,0});
        Primitive P2C(turboNormalization(0.0993707,0,1,0)*0.39195739,0,1,0,0.0993707,{0,0,0});

        Primitive P3A(turboNormalization(1.3148331,0,0,1)*0.15591627,0,0,1,1.3148331,{0,0,0});
        Primitive P3B(turboNormalization(0.3055389,0,0,1)*0.60768372,0,0,1,0.3055389,{0,0,0});
        Primitive P3C(turboNormalization(0.0993707,0,0,1)*0.39195739,0,0,1,0.0993707,{0,0,0});

        Primitive S1[3] = {S1A,S1B,S1C};
        Primitive S2[3] = {S2A,S2B,S2C};
        Primitive P1[3] = {P1A,P1B,P1C};
        Primitive P2[3] = {P2A,P2B,P2C};
        Primitive P3[3] = {P3A,P3B,P3C};

        contracted C1 (3,S1);
        contracted C2(3,S2);
        contracted C3(3,P1);
        contracted C4(3,P2);
        contracted C5(3,P3);

        basisSts.push_back(C1);
        basisSts.push_back(C2);
        basisSts.push_back(C3);
        basisSts.push_back(C4);
        basisSts.push_back(C5);
    }

    if(configuration == "He"){
        //basisSet[1];
        Nstates = 1;
        set_size(Nstates);
        Nprimitives = 3;

        //Primitive S1A(0.15432897,0,0,0,6.36242139,{0,0,0});
        //Primitive S1B(0.53532814,0,0,0,1.15892300,{0,0,0});
        //Primitive S1C(0.44463454,0,0,0,0.31364979,{0,0,0});
        //cout << turboNormalization(6.36242139,0,0,0)*0.15432897 << endl;
        Primitive S1A(0.44063,0,0,0,6.36242139,{0,0,0});
        Primitive S1B(0.426158,0,0,0,1.15892300,{0,0,0});
        Primitive S1C(0.132815,0,0,0,0.31364979,{0,0,0});

        Primitive S1[3] = {S1A,S1B,S1C};

        contracted C1 (3,S1);

        basisSts.push_back(C1);

    }
    if(configuration == "O"){
        //basisSet[3];
        Nstates = 3;
        set_size(Nstates);
        Nprimitives = 3;
        Primitive S1A(turboNormalization(130.7093200,0,0,0)*0.15432897,0,0,0,130.7093200,{0,0,0});
        Primitive S1B(turboNormalization(23.8088610,0,0,0) *0.53532814,0,0,0,23.8088610, {0,0,0});
        Primitive S1C(turboNormalization(6.4436083,0,0,0)  *0.44463454,0,0,0,6.4436083, {0,0,0});

        Primitive S2A(turboNormalization(5.0331513,0,0,0) *-0.09996723,0,0,0,5.0331513,{0,0,0});
        Primitive S2B(turboNormalization(1.1695961,0,0,0)  *0.39951283,0,0,0,1.1695961, {0,0,0});
        Primitive S2C(turboNormalization(0.3803890,0,0,0)  *0.70011547,0,0,0,0.3803890, {0,0,0});

        Primitive P1A(turboNormalization(5.0331513,1,0,0)  *0.15591627,1,0,0,5.0331513, {0,0,0});
        Primitive P1B(turboNormalization(1.1695961,1,0,0)  *0.60768372,1,0,0,1.1695961, {0,0,0});
        Primitive P1C(turboNormalization(0.3803890,1,0,0)  *0.39195739,1,0,0,0.3803890, {0,0,0});

        Primitive S1[3] = {S1A,S1B,S1C};
        Primitive S2[3] = {S2A,S2B,S2C};
        Primitive P1[3] = {P1A,P1B,P1C};

        contracted C1 (3,S1);
        contracted C2(3,S2);
        contracted C3(3,P1);

        basisSts.push_back(C1);
        basisSts.push_back(C2);
        basisSts.push_back(C3);
    }
}

double basis::nnInteraction(){
    //Contribution to the energy from the nucleon-nucleon interactions
    double result = 0;
    vec3 Rnn;
    double r;
    if(nucleusCharges.size() != 0){
        for(int i=0; i<nucleusCharges.size()-1; i++){
            for(int j=i+1; j<nucleusCharges.size(); j++){
                Rnn = nucleusPositions(i)-nucleusPositions(j);
                r = sqrt(Rnn(0)*Rnn(0)+Rnn(1)*Rnn(1)+Rnn(2)*Rnn(2));
                result += nucleusCharges(i)*nucleusCharges(j)/r;

            }
        }
    }
    return result;
}

void basis::init_integrals(){
    //Set up and solve all intergals for the current gaussian basis

    BoysFunction boys(3);
    //initializing 4 dummy primitives used in the iteration below
    Primitive A(0,0,0,0,0,{0,0,0});
    Primitive B(0,0,0,0,0,{0,0,0});
    Primitive C(0,0,0,0,0,{0,0,0});
    Primitive D(0,0,0,0,0,{0,0,0});



    for(int p=0; p<Nstates; p++){
        for(int q=0; q<Nstates; q++){
            for(int i=0; i<nPrimitivesInState(p);i++){
                for(int j=0; j<nPrimitivesInState(q);j++){
                    A = basisSts[p].getPrimitive(i);
                    B = basisSts[q].getPrimitive(j);
                    integrator AB (A,B, boys);
                    S(p,q) += AB.overlap();
                    h(p,q) += AB.kinetic();
                    for(int n = 0; n < nucleusCharges.size(); n++){
                        //add relevant interaction for each nucleus
                        AB.setupRtuv(nucleusPositions(n));
                        //cout << nucleusCharges(n)*AB.pNuclei() << endl;
                        h(p,q) -= nucleusCharges(n)*AB.pNuclei();
                    }

                    for(int r=0; r<Nstates; r++){
                        for(int s=0; s<Nstates; s++){
                            for(int k=0;k<nPrimitivesInState(r);k++){
                                for(int l=0;l<nPrimitivesInState(s);l++){

                                    C = basisSts[r].getPrimitive(k);
                                    D = basisSts[s].getPrimitive(l);
                                    //v(p,q)(r,s) += AB.pp(C,D);
                                    //v(p,r)(q,s) += AB.pp(C,D); //Changed 26/5 2014
                                    v(p,q)(r,s) += AB.pp(C,D); //Changed 27/5 2014
                                    //cout << v(p,r)(q,s) << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

double basis::factorial(double n){
    double result = 1.0;
    for(double i = 1; i< n +1; i++){
        result *= i;
    }
    return result;
}

double basis::turboNormalization(double x, double i, double j, double k){
    return pow((2*x/pi),0.75)*sqrt(pow(8*x,i+j+k)*factorial(i)*factorial(j)*factorial(k)/(factorial(2*i)*factorial(2*j)*factorial(2*k)));
}

double basis::evaluateContracted(int n, vec3 r){
    double result = 0;
    //double R = r(0)*r(0) + r(1)*r(1) + r(2) * r(2);
    double dist = 0;
    Primitive A(0,0,0,0,0,{0,0,0});
    for (int i = 0; i < Nprimitives; i++) {
        A = basisSts[n].getPrimitive(i);
        //vec3 Apos = A.nucleusPosition(1);
        vec3 dist = r-A.nucleusPosition(); //distance from nucleic center
        double dist2 = dist(0)*dist(0) + dist(1)*dist(1) + dist(2) * dist(2);
        result += A.weight()*pow(dist(0), A.xExponent())*pow(dist(1), A.yExponent())*pow(dist(2), A.zExponent())*exp(-A.exponent()*dist2);
    }
    return result;
}

void basis::set_size(int N){
    Nstates = N;
    //initializing fields and matrices
    //two-body interaction
    field<mat> vn(Nstates, Nstates);
    for (int i = 0; i < Nstates; ++i) {
        for (int j = 0; j < Nstates; ++j) {
            vn(i,j) = zeros<mat>(Nstates,Nstates); // fill V with 3x3 mx elements
        }
    }
    v = vn;

    //one-body contribution
    h.set_size(Nstates,Nstates);
    h.zeros();

    //nuclear potential
    nuclearPotential.set_size(Nstates,Nstates);
    nuclearPotential.zeros();

    //overlap matrix
    S.set_size(Nstates,Nstates);
    S.zeros();

    nPrimitivesInState.set_size(0);


}

void basis::reset(){
    for(int i = 0; i<Nstates;i++){
        nPrimitivesInState(i) = 0;
    }
    nPrimitivesInState.set_size(0);
    nPrimitivesInState = vec ();
    basisSts.resize(3);
    set_size(0);
}

void basis::printAllContracted(){
    Primitive A(0,0,0,0,0,{0,0,0});
    for(int n=0;n<Nstates;n++){
        for(int i = 0; i<Nprimitives;i++){
            A = basisSts[n].getPrimitive(i);
            cout << "Contracted/primitive:" << n << "/" << i << ":" << A.exponent() << " " << A.weight() << " " << A.xExponent() << A.yExponent() << A.zExponent() << endl;
        }
        cout << endl;
    }
}

void basis::init_HTO4(double nProtons){
    //initialize 4 hydrogenlike orbitals
    nucleusPositions.set_size(1);
    nucleusCharges.set_size(1);
    nucleusPositions(0) = {0,0,0};
    nucleusCharges(0) = (double) nProtons;

    set_size(4);
    set_orthonormal();
    Z = (double) nProtons;
    for(int i = 0; i<Nstates; i++){
        h(i,i) = h0(i,i);
    }

    v(0,0)(0,0)= 5*Z/8;
    v(0,0)(0,1)= 4096*sqrt(2)*Z/64827;
    v(0,0)(0,2)= 1269*sqrt(3)*Z/50000;
    v(0,0)(0,3)= 416415744*Z/15083778125;
    v(0,0)(1,0)= 4096*sqrt(2)*Z/64827;
    v(0,0)(1,1)= 16*Z/729;
    v(0,0)(1,2)= 110592*sqrt(6)*Z/24137569;
    v(0,0)(1,3)= 98304*sqrt(2)*Z/19487171;
    v(0,0)(2,0)= 1269*sqrt(3)*Z/50000;
    v(0,0)(2,1)= 110592*sqrt(6)*Z/24137569;
    v(0,0)(2,2)= 189*Z/32768;
    v(0,0)(2,3)= 1808400384*sqrt(3)*Z/852891037441;
    v(0,0)(3,0)= 416415744*Z/15083778125;
    v(0,0)(3,1 )= 98304*sqrt(2)*Z/19487171;
    v(0,0)(3,2 )= 1808400384*sqrt(3)*Z/852891037441;
    v(0,0)( 3,3)= 22848*Z/9765625;
    v(0,1)( 0,0 )= 4096*sqrt(2)*Z/64827;
    v(0,1)( 0,1 )= 17*Z/81;
    v(0,1)( 0,2 )= 1555918848*sqrt(6)*Z/75429903125;
    v(0,1)( 0,3 )= 288456704*sqrt(2)*Z/14206147659;
    v(0,1)( 1,0 )= 16*Z/729;
    v(0,1)( 1,1 )= 512*sqrt(2)*Z/84375;
    v(0,1)( 1,2 )= 2160*sqrt(3)*Z/823543;
    v(0,1)( 1,3 )= 376832*Z/129140163;
    v(0,1)( 2,0)= 110592*sqrt(6)*Z/24137569;
    v(0,1)( 2,1 )= 29943*sqrt(3)*Z/13176688;
    v(0,1)( 2,2 )= 1216512*sqrt(2)*Z/815730721;
    v(0,1)( 2,3 )= 423788544*sqrt(6)*Z/762939453125;
    v(0,1)( 3,0 )= 98304*sqrt(2)*Z/19487171;
    v(0,1)( 3,1 )= 108838912*Z/44840334375;
    v(0,1)( 3,2 )= 10148806656*sqrt(6)*Z/19073486328125;
    v(0,1)(3,3)= 39*Z/32768*sqrt(2);
    v(0,2)(0,0 )= 1269*sqrt(3)*Z/50000;
    v(0,2)(0,1 )= 1555918848*sqrt(6)*Z/75429903125;
    v(0,2)( 0,2 )= 815*Z/8192;
    v(0,2)( 0,3 )= 11694850770862080*sqrt(3)*Z/702392443647273463;
    v(0,2)( 1,0 )= 110592*sqrt(6)*Z/24137569;
    v(0,2)( 1,1 )= 2160*sqrt(3)*Z/823543;
    v(0,2)( 1,2 )= 37826560*sqrt(2)*Z/22024729467;
    v(0,2)( 1,3 )= 487489536*sqrt(6)*Z/762939453125;
    v(0,2)( 2,0)= 189*Z/32768;
    v(0,2)( 2,1 )= 1216512*sqrt(2)*Z/815730721;
    v(0,2)( 2,2 )= 617*Z/314928*sqrt(3);
    v(0,2)( 2,3 )= 30254432256*Z/41426511213649;
    v(0,2)( 3,0 )= 1808400384*sqrt(3)*Z/852891037441;
    v(0,2)( 3,1 )= 10148806656*sqrt(6)*Z/19073486328125;
    v(0,2)( 3,2 )= 90581886173184*Z/129457847542653125;
    v(0,2)( 3,3)= 74450880*sqrt(3)*Z/285311670611;
    v(0,3)( 0,0 )= 416415744*Z/15083778125;
    v(0,3)( 0,1 )= 288456704*sqrt(2)*Z/14206147659;
    v(0,3)( 0,2 )= 11694850770862080*sqrt(3)*Z/702392443647273463;
    v(0,3)( 0,3 )= 22513*Z/390625;
    v(0,3)( 1,0 )= 98304*sqrt(2)*Z/19487171;
    v(0,3)( 1,1 )= 376832*Z/129140163;
    v(0,3)( 1,2 )= 487489536*sqrt(6)*Z/762939453125;
    v(0,3)( 1,3 )= 5053*Z/3538944*sqrt(2);
    v(0,3)( 2,0)= 1808400384*sqrt(3)*Z/852891037441;
    v(0,3)( 2,1 )= 423788544*sqrt(6)*Z/762939453125;
    v(0,3)( 2,2 )= 30254432256*Z/41426511213649;
    v(0,3)( 2,3 )= 1243165779*sqrt(3)*Z/4564986729776;
    v(0,3)( 3,0 )= 22848*Z/9765625;
    v(0,3)( 3,1 )= 39*Z/32768*sqrt(2);
    v(0,3)( 3,2 )= 74450880*sqrt(3)*Z/285311670611;
    v(0,3)( 3,3)= 1804351488*Z/6179146071875;
    v(1,0)( 0,0 )= 4096*sqrt(2)*Z/64827;
    v(1,0)( 0,1 )= 16*Z/729;
    v(1,0)( 0,2 )= 110592*sqrt(6)*Z/24137569;
    v(1,0)( 0,3 )= 98304*sqrt(2)*Z/19487171;
    v(1,0)(1,0 )= 17*Z/81;
    v(1,0)(1,1 )= 512*sqrt(2)*Z/84375;
    v(1,0)( 1,2 )= 29943*sqrt(3)*Z/13176688;
    v(1,0)( 1,3 )= 108838912*Z/44840334375;
    v(1,0)( 2,0)= 1555918848*sqrt(6)*Z/75429903125;
    v(1,0)( 2,1 )= 2160*sqrt(3)*Z/823543;
    v(1,0)( 2,2 )= 1216512*sqrt(2)*Z/815730721;
    v(1,0)( 2,3 )= 10148806656*sqrt(6)*Z/19073486328125;
    v(1,0)(3,0 )= 288456704*sqrt(2)*Z/14206147659;
    v(1,0)(3,1 )= 376832*Z/129140163;
    v(1,0)(3,2 )= 423788544*sqrt(6)*Z/762939453125;
    v(1,0)(3,3)= 39*Z/32768*sqrt(2);
    v(1,1)(0,0 )= 16*Z/729;
    v(1,1)(0,1 )= 512*sqrt(2)*Z/84375;
    v(1,1)(0,2 )= 2160*sqrt(3)*Z/823543;
    v(1,1)(0,3 )= 376832*Z/129140163;
    v(1,1)(1,0 )= 512*sqrt(2)*Z/84375;
    v(1,1)(1,1 )= 77*Z/512;
    v(1,1)(1,2 )= 5870679552*sqrt(6)*Z/669871503125;
    v(1,1)(1,3 )= 31363072*sqrt(2)*Z/4202539929;
    v(1,1)(2,0)= 2160*sqrt(3)*Z/823543;
    v(1,1)(2,1 )= 5870679552*sqrt(6)*Z/669871503125;
    v(1,1)(2,2 )= 73008*Z/9765625;
    v(1,1)(2,3 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(1,1)( 3,0 )= 376832*Z/129140163;
    v(1,1)( 3,1 )= 31363072*sqrt(2)*Z/4202539929;
    v(1,1)( 3,2 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(1,1)( 3,3)= 424*Z/177147;
    v(1,2)( 0,0 )= 110592*sqrt(6)*Z/24137569;
    v(1,2)( 0,1 )= 2160*sqrt(3)*Z/823543;
    v(1,2)( 0,2 )= 37826560*sqrt(2)*Z/22024729467;
    v(1,2)( 0,3 )= 487489536*sqrt(6)*Z/762939453125;
    v(1,2)( 1,0 )= 29943*sqrt(3)*Z/13176688;
    v(1,2)( 1,1 )= 5870679552*sqrt(6)*Z/669871503125;
    v(1,2)( 1,2 )= 32857*Z/390625;
    v(1,2)( 1,3 )= 55508689880137728*sqrt(3)*Z/5049196699148208943;
    v(1,2)( 2,0)= 1216512*sqrt(2)*Z/815730721;
    v(1,2)( 2,1 )= 73008*Z/9765625;
    v(1,2)( 2,2 )= 6890942464*sqrt(2./3)*Z/1210689028125;
    v(1,2)( 2,3 )= 69158928384*sqrt(2.)*Z/34271896307633;
    v(1,2)( 3,0 )= 423788544*sqrt(6)*Z/762939453125;
    v(1,2)( 3,1 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(1,2)( 3,2 )= 36645380390912*sqrt(2)*Z/24984212408264457;
    v(1,2)( 3,3)= 145503*sqrt(3/2)*Z/134217728;
    v(1,3)( 0,0 )= 98304*sqrt(2)*Z/19487171;
    v(1,3)( 0,1 )= 376832*Z/129140163;
    v(1,3)(0,2 )= 487489536*sqrt(6)*Z/762939453125;
    v(1,3)( 0,3 )= 5053*Z/3538944*sqrt(2);
    v(1,3)( 1,0 )= 108838912*Z/44840334375;
    v(1,3)( 1,1 )= 31363072*sqrt(2)*Z/4202539929;
    v(1,3)( 1,2 )= 55508689880137728*sqrt(3)*Z/5049196699148208943;
    v(1,3)( 1,3 )= 4043*Z/78732;
    v(1,3)( 2,0)= 10148806656*sqrt(6)*Z/19073486328125;
    v(1,3)( 2,1 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(1,3)( 2,2 )= 69158928384*sqrt(2)*Z/34271896307633;
    v(1,3)( 2,3 )= 2496169683*sqrt(3/2)*Z/1677721600000;
    v(1,3)( 3,0 )= 39*Z/32768*sqrt(2);
    v(1,3)( 3,1 )= 424*Z/177147;;
    v(1,3)( 3,2 )= 145503*sqrt(3/2)*Z/134217728;;
    v(1,3)( 3,3)= 21252608*sqrt(2)*Z/35595703125;
    v(2,0)( 0,0 )= 1269*sqrt(3)*Z/50000;
    v(2,0)( 0,1 )= 110592*sqrt(6)*Z/24137569;
    v(2,0)( 0,2 )= 189*Z/32768;
    v(2,0)( 0,3 )= 1808400384*sqrt(3)*Z/852891037441;
    v(2,0)( 1,0 )= 1555918848*sqrt(6)*Z/75429903125;
    v(2,0)( 1,1 )= 2160*sqrt(3)*Z/823543;
    v(2,0)( 1,2 )= 1216512*sqrt(2)*Z/815730721;
    v(2,0)( 1,3 )= 10148806656*sqrt(6)*Z/19073486328125;
    v(2,0)( 2,0)= 815*Z/8192;
    v(2,0)( 2,1 )= 37826560*sqrt(2)*Z/22024729467;
    v(2,0)( 2,2 )= 617*Z/314928*sqrt(3);
    v(2,0)( 2,3 )= 90581886173184*Z/129457847542653125;
    v(2,0)( 3,0 )= 11694850770862080*sqrt(3)*Z/702392443647273463;
    v(2,0)( 3,1 )= 487489536*sqrt(6)*Z/762939453125;
    v(2,0)( 3,2 )= 30254432256*Z/41426511213649;
    v(2,0)( 3,3)= 74450880*sqrt(3)*Z/285311670611;
    v(2,1)( 0,0 )= 110592*sqrt(6)*Z/24137569;
    v(2,1)( 0,1 )= 29943*sqrt(3)*Z/13176688;
    v(2,1)( 0,2 )= 1216512*sqrt(2)*Z/815730721;
    v(2,1)( 0,3 )= 423788544*sqrt(6)*Z/762939453125;
    v(2,1)( 1,0 )= 2160*sqrt(3)*Z/823543;
    v(2,1)( 1,1 )= 5870679552*sqrt(6)*Z/669871503125;
    v(2,1)( 1,2 )= 73008*Z/9765625;
    v(2,1)( 1,3 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(2,1)( 2,0)= 37826560*sqrt(2)*Z/22024729467;
    v(2,1)(2,1 )= 32857*Z/390625;
    v(2,1)(2,2 )= 6890942464*sqrt(2./3.)*Z/1210689028125;
    v(2,1)(2,3 )= 36645380390912*sqrt(2)*Z/24984212408264457;
    v(2,1)(3,0 )= 487489536*sqrt(6)*Z/762939453125;
    v(2,1)(3,1 )= 55508689880137728*sqrt(3)*Z/5049196699148208943;
    v(2,1)(3,2 )= 69158928384*sqrt(2)*Z/34271896307633;
    v(2,1)( 3,3)= 145503*sqrt(3/2)*Z/134217728;
    v(2,2)( 0,0 )= 189*Z/32768;
    v(2,2)(0,1 )= 1216512*sqrt(2)*Z/815730721;
    v(2,2)( 0,2 )= 617*Z/314928*sqrt(3);
    v(2,2)( 0,3 )= 30254432256*Z/41426511213649;
    v(2,2)( 1,0 )= 1216512*sqrt(2)*Z/815730721;
    v(2,2)( 1,1 )= 73008*Z/9765625;
    v(2,2)( 1,2 )= 6890942464*sqrt(2./3)*Z/1210689028125;
    v(2,2)(1,3 )= 69158928384*sqrt(2)*Z/34271896307633;
    v(2,2)( 2,0)= 617*Z/314928*sqrt(3);
    v(2,2)( 2,1 )= 6890942464*sqrt(2./3)*Z/1210689028125;
    v(2,2)( 2,2 )= 17*Z/256;
    v(2,2)( 2,3 )= 2486755845603328*Z/158298797548828125*sqrt(3);
    v(2,2)( 3,0 )= 30254432256*Z/41426511213649;
    v(2,2)( 3,1 )= 69158928384*sqrt(2)*Z/34271896307633;
    v(2,2)( 3,2 )= 2486755845603328*Z/158298797548828125*sqrt(3);
    v(2,2)( 3,3)= 2560158144*Z/678223072849;
    v(2,3)( 0,0 )= 1808400384*sqrt(3)*Z/852891037441;
    v(2,3)( 0,1 )= 423788544*sqrt(6)*Z/762939453125;
    v(2,3)( 0,2 )= 30254432256*Z/41426511213649;
    v(2,3)( 0,3 )= 1243165779*sqrt(3)*Z/4564986729776;
    v(2,3)( 1,0 )= 10148806656*sqrt(6)*Z/19073486328125;
    v(2,3)( 1,1 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(2,3)( 1,2 )= 69158928384*sqrt(2)*Z/34271896307633;
    v(2,3)( 1,3 )= 2496169683*sqrt(3/2)*Z/1677721600000;
    v(2,3)(2,0)= 90581886173184*Z/129457847542653125;
    v(2,3)( 2,1 )= 36645380390912*sqrt(2)*Z/24984212408264457;
    v(2,3)( 2,2 )= 2486755845603328*Z/158298797548828125*sqrt(3);
    v(2,3)( 2,3 )= 621550729*Z/13841287201;
    v(2,3)( 3,0 )= 74450880*sqrt(3)*Z/285311670611;
    v(2,3)( 3,1 )= 145503*sqrt(3/2)*Z/134217728;
    v(2,3)( 3,2 )= 2560158144*Z/678223072849;
    v(2,3)(3,3)= 413631006610176000.0*sqrt(3.0)*Z/249430673908303812379.0;
    v(3,0)( 0,0 )= 416415744*Z/15083778125;
    v(3,0)( 0,1 )= 98304*sqrt(2)*Z/19487171;
    v(3,0)( 0,2 )= 1808400384*sqrt(3)*Z/852891037441;
    v(3,0)( 0,3 )= 22848*Z/9765625;
    v(3,0)( 1,0 )= 288456704*sqrt(2)*Z/14206147659;
    v(3,0)( 1,1 )= 376832*Z/129140163;
    v(3,0)( 1,2 )= 423788544*sqrt(6)*Z/762939453125;
    v(3,0)( 1,3 )= 39*Z/32768*sqrt(2);
    v(3,0)( 2,0)= 11694850770862080*sqrt(3)*Z/702392443647273463;
    v(3,0)( 2,1 )= 487489536*sqrt(6)*Z/762939453125;
    v(3,0)( 2,2 )= 30254432256*Z/41426511213649;
    v(3,0)( 2,3 )= 74450880*sqrt(3)*Z/285311670611;
    v(3,0)( 3,0 )= 22513*Z/390625;
    v(3,0)( 3,1 )= 5053*Z/3538944*sqrt(2);
    v(3,0)( 3,2 )= 1243165779*sqrt(3)*Z/4564986729776;
    v(3,0)( 3,3)= 1804351488*Z/6179146071875;
    v(3,1)( 0,0 )= 98304*sqrt(2)*Z/19487171;
    v(3,1)( 0,1 )= 108838912*Z/44840334375;
    v(3,1)( 0,2 )= 10148806656*sqrt(6)*Z/19073486328125;
    v(3,1)( 0,3 )= 39*Z/32768*sqrt(2);
    v(3,1)( 1,0 )= 376832*Z/129140163;
    v(3,1)( 1,1 )= 31363072*sqrt(2)*Z/4202539929;
    v(3,1)( 1,2 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(3,1)( 1,3 )= 424*Z/177147;
    v(3,1)( 2,0)= 487489536*sqrt(6)*Z/762939453125;
    v(3,1)( 2,1 )= 55508689880137728*sqrt(3)*Z/5049196699148208943;
    v(3,1)(2,2 )= 69158928384*sqrt(2)*Z/34271896307633;
    v(3,1)( 2,3 )= 145503*sqrt(3/2)*Z/134217728;
    v(3,1)( 3,0 )= 5053*Z/3538944*sqrt(2);
    v(3,1)( 3,1 )= 4043*Z/78732;
    v(3,1)( 3,2 )= 2496169683*sqrt(3/2)*Z/1677721600000;
    v(3,1)( 3,3)= 21252608*sqrt(2)*Z/35595703125;
    v(3,2)(0,0 )= 1808400384*sqrt(3)*Z/852891037441;
    v(3,2)( 0,1 )= 10148806656*sqrt(6)*Z/19073486328125;
    v(3,2)( 0,2 )= 90581886173184*Z/129457847542653125;
    v(3,2)( 0,3 )= 74450880*sqrt(3)*Z/285311670611;
    v(3,2)( 1,0 )= 423788544*sqrt(6)*Z/762939453125;
    v(3,2)( 1,1 )= 14739259392*sqrt(3)*Z/6131066257801;
    v(3,2)(1,2 )= 36645380390912*sqrt(2)*Z/24984212408264457;
    v(3,2)( 1,3 )= 145503*sqrt(3/2)*Z/134217728;
    v(3,2)( 2,0)= 30254432256*Z/41426511213649;
    v(3,2)( 2,1 )= 69158928384*sqrt(2)*Z/34271896307633;
    v(3,2)( 2,2 )= 2486755845603328*Z/158298797548828125*sqrt(3);
    v(3,2)( 2,3 )= 2560158144*Z/678223072849;
    v(3,2)( 3,0 )= 1243165779*sqrt(3)*Z/4564986729776;
    v(3,2)( 3,1 )= 2496169683*sqrt(3/2)*Z/1677721600000;
    v(3,2)( 3,2 )= 621550729*Z/13841287201;
    v(3,2)( 3,3)= 413631006610176000.0*sqrt(3.0)*Z/249430673908303812379.0;
    v(3,3)( 0,0 )= 22848*Z/9765625;
    v(3,3)( 0,1 )= 39*Z/32768*sqrt(2);
    v(3,3)( 0,2 )= 74450880*sqrt(3)*Z/285311670611;
    v(3,3)( 0,3 )= 1804351488*Z/6179146071875;
    v(3,3)( 1,0 )= 39*Z/32768*sqrt(2);
    v(3,3)( 1,1 )= 424*Z/177147;
    v(3,3)( 1,2 )= 145503*sqrt(3/2)*Z/134217728;
    v(3,3 )(1,3 )= 21252608*sqrt(2)*Z/35595703125;
    v(3,3)( 2,0)= 74450880*sqrt(3)*Z/285311670611;
    v(3,3)( 2,1 )= 145503*sqrt(3/2)*Z/134217728;
    v(3,3)( 2,2 )= 2560158144*Z/678223072849;
    v(3,3)( 2,3 )= 413631006610176000.0*sqrt(3.0)*Z/249430673908303812379.0;
    v(3,3)( 3,0 )= 1804351488*Z/6179146071875;
    v(3,3)( 3,1 )= 21252608*sqrt(2)*Z/35595703125;
    v(3,3)( 3,2 )= 413631006610176000*sqrt(3.0)*Z/249430673908303812379.0;
    v(3,3)(3,3)= 19541*Z/524288;


    flip_index();
}

void basis::flip_index(){
    field<mat> vn(Nstates, Nstates);
    for (int p = 0; p < Nstates; ++p) {
        for (int q = 0; q < Nstates; ++q) {
            vn(p,q) = zeros<mat>(Nstates,Nstates); // fill V with 3x3 mx elements
        }
    }
    for (int p = 0; p < Nstates; ++p) {
        for (int q = 0; q < Nstates; ++q) {
            for(int r=0;r<Nstates;++r){
                for(int s=0;s<Nstates;++s){
                    vn(p,r)(q,s) = v(p,q)(r,s);
                }
            }
        }
    }
    v = vn;
}

