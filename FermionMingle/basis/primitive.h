#ifndef PRIMITIVE_H
#define PRIMITIVE_H
#include <armadillo>
using namespace arma;

class Primitive {

public:
    explicit Primitive(double weight, int xExponent,
                       int yExponent, int zExponent, double exponent,
                       vec nucleusPosition);
/*
    inline double exponent() const;
    inline int zExponent() const;
    inline int yExponent() const;
    inline int xExponent() const;
    inline double weight() const;
    inline const vec& nucleusPosition() const;
*/
    double exponent() const;
    int zExponent() const;
    int yExponent() const;
    int xExponent() const;
    double weight() const;
    const vec& nucleusPosition() const;

private:
    double m_weight;
    int m_xExponent;
    int m_yExponent;
    int m_zExponent;
    double m_exponent;
    vec m_nucleusPosition;
};

// this is a way much faster implementation of the return funcitons, using inline
// The functions it then "glued" directly into the code where they are called,
// this way, we skip an extra call for each function.
// the back side of this is if something goes wrong. It can be painfull to find the
// error ! NOTE: we do not declare the inline funcitons in the .cpp file!

/*  Therefore we should test out if the class works the usual way first!

inline double Primitive::exponent() const
{
return m_exponent;
}
inline int Primitive::zExponent() const
{
return m_zExponent;
}
inline int Primitive::yExponent() const
{
inline return m_yExponent;
}
inline int Primitive::xExponent() const
{
inline return m_xExponent;
}
inline double Primitive::weight() const
{
inline return m_weight;
}
inline const vec& Primitive::nucleusPosition() const
{
return m_nucleusPosition;
}
*/

#endif // PRIMITIVE_H
