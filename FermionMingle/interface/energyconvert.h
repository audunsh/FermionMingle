#ifndef ENERGYCONVERT_H
#define ENERGYCONVERT_H

#include <string>

using namespace std;

class energyconvert
{
public:
    energyconvert();
    energyconvert(double E0);
    double as(string units);
    double E;
};

#endif // ENERGYCONVERT_H
