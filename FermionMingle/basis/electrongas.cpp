#include "electrongas.h"
#include <armadillo>


using namespace std;
using namespace arma;

electrongas::electrongas()
{
}

void electrongas::generate_state_list(int Ne, double rs){
    N = Ne;

    //Volum = nokkuperte*4.d0*pi*r_s**3/3.d0
    r_s = rs;


    //prefactor1 = 4*pi/(L*L*L); //These are not necessarily correct

    prefactor1 = 3/(14.0*r_s*r_s*r_s);
    prefactor2 = .5;

    double Nmax = sqrt(N) + 1;
    int energy = 0;
    int nStates = 0;

    //Counting the number of states needed up to energy level N
    for(int x = -Nmax; x<Nmax+1; x++){
        for(int y = -Nmax; y<Nmax+1; y++){
            for(int z = -Nmax; z<Nmax+1; z++){
                energy = x*x + y*y + z*z;
                if(energy < N + 1){
                    //cout << "Energy:" << energy << " State: " << x << " " << y << " " <<  z << endl;
                    nStates += 2; //Due to spin degeneracy
                }
            }
        }
    }
    //cout << "Up to energy level " << N << " there will be " << nStates << " states." << endl;
    prefactor1 = 3.0/(nStates*r_s*r_s*r_s);
    L3 = nStates*4.0*pi*r_s*r_s*r_s/3.0;
    L = pow(L3, 1.0/3.0);
    //Setting up all all states
    mat k_combinations = zeros(nStates, 5);
    sorted_energy.zeros(nStates, 4);
    int index_count = 0;
    double e2;

    for(int x = -Nmax; x<Nmax+1; x++){
        for(int y = -Nmax; y<Nmax+1; y++){
            for(int z = -Nmax; z<Nmax+1; z++){
                energy = (x*x + y*y + z*z);
                e2 = energy*pi*pi/(2*L*L);
                if(energy < N + 1){
                    k_combinations(index_count, 0) = e2; //energy*prefactor2*(53.63609*pi*pi/L3);
                    k_combinations(index_count, 1) = x;
                    k_combinations(index_count, 2) = y;
                    k_combinations(index_count, 3) = z;
                    k_combinations(index_count, 4) = 0;
                    index_count += 1;

                    k_combinations(index_count, 0) = e2; //energy*prefactor2*(53.63609*pi*pi/L3);
                    k_combinations(index_count, 1) = x;
                    k_combinations(index_count, 2) = y;
                    k_combinations(index_count, 3) = z;
                    k_combinations(index_count, 4) = -1;
                    index_count += 1;
                }
            }
        }
    }
    //k_combinations.print();

    vec temp_vec = k_combinations.col(0);
    Energy.zeros(index_count);
    n_basis_functions = index_count;
    uvec sorted_vector = sort_index(temp_vec);

    for(int i = 0; i< index_count; i++){
        for(int j = 0; j< 4; j++){
            sorted_energy(i,j) = k_combinations(sorted_vector(i), j+1);
        }
        Energy(i) = k_combinations(sorted_vector(i), 0);
    }
    //sorted_energy.print();
    Energy.print();
}

double electrongas::absdiff2(vec A, vec B){
    double D = 0;
    for (int i =0; i < 3; i++){
        D += (A(i) - B(i))*(A(i) - B(i));
    }
    return D;
}

int electrongas::kd(int A, int B){
    return 1*(A==B);
}

int electrongas::kd_vec(rowvec A, rowvec B){
    int D = 1;
    for(int i = 0; i < A.n_elem; i++){
        D*=(A(i)==B(i));
    }
    return D;
}

double electrongas::h(int P, int Q){
    return Energy(P)*(P==Q);
}

double electrongas::v(int P, int Q, int R, int S){
    //Two body interaction
    rowvec KP = sorted_energy.row(P);
    rowvec KQ = sorted_energy.row(Q);
    rowvec KR = sorted_energy.row(R);
    rowvec KS = sorted_energy.row(S);

    //Two electron interaction
    double value;
    double term1= 0;
    double term2= 0;
    double kd1, kd2;
    int spinP = KP(3);
    int spinQ = KQ(3);
    int spinR = KR(3);
    int spinS = KS(3);

    rowvec kp, kq, kr, ks;
    kp << KP(0) << KP(1) << KP(2);
    kq << KQ(0) << KQ(1) << KQ(2);
    kr << KR(0) << KR(1) << KR(2);
    ks << KS(0) << KS(1) << KS(2);

    value = kd_vec((kp+kq), (kr+ks));
    if(value ==0 ){
        return 0;
    }
    else{
        term1 = kd(spinP, spinR)*kd(spinQ,spinS);
        kd1 = kd_vec(kp, kr);
        if(kd1!= 1.0){
            term1 = term1 / absdiff2(kr, kp);
        }

        term2 = kd(spinP, spinS)*kd(spinQ,spinR);
        kd2 = kd_vec(kp, ks);
        if(kd2!= 1.0){
            term2 = term2 / absdiff2(ks, kp);
        }

        value *= (term1 - term2);
        return prefactor1*value;
    }
}

double electrongas::f(int P, int Q){
    //Fock operator matrix elements
    double val = prefactor2;
    rowvec KP = sorted_energy.row(P);
    rowvec KQ = sorted_energy.row(Q);
    vec kp;
    kp.zeros(3);
    kp(0) = KP(0);
    kp(1) = KP(1);
    kp(2) = KP(2);
    val *= dot(kp, kp);
    val *= kd_vec(KP,KQ);
    double val2 = 0;
    for(int i = 0; i < n_basis_functions; i++){
        if((i != P) && (i != Q)){
            val2 += v(P, i, Q, i);
        }
    }
    val += prefactor1*val2;
    return val;
}
