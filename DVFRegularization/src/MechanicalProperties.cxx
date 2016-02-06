#include "MechanicalProperties.h"

using namespace dvfRegularization;

MechanicalProperties::MechanicalProperties(float youngsModulus, float poissonRatio) 
: m_YoungsModulus(youngsModulus), m_PoissonRatio(poissonRatio)
{
}

MechanicalProperties::~MechanicalProperties(void)
{
}
