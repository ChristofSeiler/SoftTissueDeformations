/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Regularize vector displacement field.
 */

#include "APrioriImageInformation.h"

APrioriImageInformation::APrioriImageInformation(float deformation, 
												 int directionOfDeformation,
												 int boundaryCondition) 
												 : m_deformation(deformation), 
												 m_directionOfDeformation(directionOfDeformation),
												 m_boundaryCondition(boundaryCondition) {
}

APrioriImageInformation::~APrioriImageInformation(void) {
}

void APrioriImageInformation::computeInformation(void) {
	// extract information out of image
}

float APrioriImageInformation::getLocalDeformationForTissueType(int tissueType) {
	return m_localTissueDeformation[tissueType];
}