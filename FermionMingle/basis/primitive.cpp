#include "basis/primitive.h"

Primitive::Primitive(double weight,
                     int xExponent, int yExponent, int zExponent,
                     double exponent, vec nucleusPosition) :

m_weight(weight),
m_xExponent(xExponent),
m_yExponent(yExponent),
m_zExponent(zExponent),
m_exponent(exponent),
m_nucleusPosition(nucleusPosition)

{

}

double Primitive::exponent() const
{
return m_exponent;
}
int Primitive::zExponent() const
{
return m_zExponent;
}
int Primitive::yExponent() const
{
return m_yExponent;
}
int Primitive::xExponent() const
{
return m_xExponent;
}
double Primitive::weight() const
{
return m_weight;
}
const vec& Primitive::nucleusPosition() const
{
return m_nucleusPosition;
}
