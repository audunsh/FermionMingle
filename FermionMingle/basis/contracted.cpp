#include "basis/contracted.h"
#include <armadillo>
#include "basis/primitive.h"

using namespace std;
using namespace arma;

contracted::contracted(){
    Nprimitives = 0;
}

contracted::contracted(int N, Primitive primitives[]){
    Nprimitives = N;
    //basisFunction[N];
    //vector basisFs.resize(N);
    //basisFs.resize(N);
    for(int i=0; i<N; i++){
        //basisFunction[i] = primitives[i];
        //basisFs.push_back(primitives[i]);
        appendPrimitive(primitives[i]);

    }
}

void contracted::appendPrimitive(Primitive P){
    basisFs.push_back(P);
}

void contracted::setPrimitive(int n){
}

Primitive contracted::getPrimitive(int n){
    //return basisFunction[n];
    return basisFs[n];
}

void contracted::free(){
    delete[] basisFunction;
}
