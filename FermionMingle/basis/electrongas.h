#ifndef ELECTRONGAS_H
#define ELECTRONGAS_H
#include <armadillo>


using namespace std;
using namespace arma;

class electrongas
{
public:
    electrongas();
    void generate_state_list(int Ne, int Le);
    double absdiff2(vec A, vec B);
    double v(int P, int Q, int R, int S);
    double f(int P, int Q);
    int kd_vec(rowvec A, rowvec B);
    int kd(int A, int B);

    int L;
    int N;
    int n_basis_functions;

    mat sorted_energy;
    vec Energy;
    double prefactor1;
    double prefactor2;

    double pi = 4*atan(1);
};

#endif // ELECTRONGAS_H
