#ifndef ENERGYCONVERT_H
#define ENERGYCONVERT_H

#include <string>

using namespace std;

class energyconvert
{
public:
    energyconvert();
    energyconvert(double E0);
    double E;
    double hartree();
    double angstrom();
    double au();
    double rydberg();
    double joule();
    double eV();
};

#endif // ENERGYCONVERT_H
