/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Regularize vector displacement field.
 */

#pragma once

#include "Definitions.h"

class APrioriImageInformation
{
public:
	APrioriImageInformation(float deformation, int directionOfDeformation, int boundaryCondition);
	~APrioriImageInformation(void);

	void computeInformation(void);
	float getLocalDeformationForTissueType(int tissueType);

	float getDeformation(void) { return m_deformation; }
	void setDeformation(float deformation) { m_deformation = deformation; }
	int getDirectionOfDeformation(void);
	void setDirectionOfDeformation(int directionOfDeformation);
	int getBoundaryCondition(void) { return m_boundaryCondition; }
	void setBoundaryCondition(int boundaryCondition) { m_boundaryCondition = boundaryCondition; }
private:
	float m_deformation;
	int m_boundaryCondition;
	int m_directionOfDeformation;
	std::vector<float> m_localTissueDeformation;
};
