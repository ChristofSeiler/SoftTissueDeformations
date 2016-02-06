/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Implementation of Metropolis Algorithm for sampling the associated Gibbs 
 *			distribution and minimizing the energy function.
 */

#pragma once

#include "Definitions.h"
#include "MechanicalProperties.h"

namespace dvfRegularization {

	class MetropolisSamplerMinimizer {
	public:
		void run(VectorImageType* inital, std::vector<ImageType::Pointer>& elasticVector, ImageType* label, 
			std::vector<MechanicalProperties*>& mechanicalProperties, VectorImageType* abaqus);
	};
}