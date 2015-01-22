#include "energyconvert.h"
#include <string>
#include <iostream>

using namespace std;

energyconvert::energyconvert()
{
}

energyconvert::energyconvert(double E0){
    E = E0; //Energy in atomic units (hartrees)
}

double energyconvert::rydberg(){return E*2;}

double energyconvert::au(){return E;}

double energyconvert::hartree(){return E;}

double energyconvert::angstrom(){return E*0.529177;}

double energyconvert::eV(){return E*27.212;}

double energyconvert::joule(){return E*4.359810E-18;}
