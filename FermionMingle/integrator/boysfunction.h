#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include <armadillo>
#include <string>
#include <fstream>

using namespace arma;
using namespace std;

class BoysFunction
{
public:
    BoysFunction();
    BoysFunction(int angMomMax);
    void setx(double x);
    double returnValue(int n);

private:
    double tabulated(int n, double x);
    double asymptotic(int n, double x);
    double factorial2(int n);

    mat m_Ftabulated;
    vec m_F;
    int m_nMax;
};

#endif // BOYSFUNCTION_H

