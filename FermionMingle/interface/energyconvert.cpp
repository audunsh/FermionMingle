#include "energyconvert.h"
#include <string>

using namespace std;

energyconvert::energyconvert()
{
}

energyconvert::energyconvert(double E0){
    E = E0; //Energy in atomic units (hartrees)
}

double energyconvert::as(string units){
    if(units == "rydberg" or "Rydberg" or "Rydbergs" or "rydbergs" or "hartree" or "hartrees" or "Hartree" or "Hartrees" or "au" or "atomic units" or "bohr" or "Bohr" or "angstrom" or "Angstrom" or "ev" or "eV" or "electron volt" or "joule" or "J" or "Joule"){
        if(units == "hartree" or "hartrees" or "Hartree" or "Hartrees" or "au" or "atomic units" or "bohr" or "Bohr"){
            return E;
        }
        if(units == "angstrom" or "Angstrom"){
            return E*0.529177;
        }
        if(units == "rydberg" or "Rydberg" or "Rydbergs" or "rydbergs"){
            return E*2;
        }
        if(units == "ev" or "eV" or "electron volt" or "Electron Volt"){
            return E*27.212;
        }
        if(units == "joule" or "J" or "Joule"){
            return E*4.359810E-18;
        }


    }
    else{
        return 0;
    }

}
