/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Mechanical properties for a certain type of pixel values.
 */

#pragma once

namespace dvfRegularization {

class MechanicalProperties
{
public:
	MechanicalProperties(float youngsModulus, float poissonRatio);
	~MechanicalProperties(void);
	void setYoungsModulus(float youngsModulus) { m_YoungsModulus = youngsModulus; }
	float getYoungsModulus() { return m_YoungsModulus; }
	void setPoissonRatio(float poissonRatio) { m_PoissonRatio = poissonRatio; }
	float getPoissonRatio() { return m_PoissonRatio; }
private:
	float m_YoungsModulus;
	float m_PoissonRatio;
};

}