/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Represents an energy function with its two parts.
 *			Umodel + Uprior = Utotal
 */

#ifndef _itkEnergyFunction_txx
#define _itkEnergyFunction_txx
#include "itkEnergyFunction.h"

#include "itkNumericTraits.h"

namespace itk {

template<class TImage,class TNLabelImage>
EnergyFunction<TImage,TNLabelImage>::EnergyFunction():
	m_Alpha(0.0), m_Beta(0.0), m_Gamma(0.0)
{
}

template<class TImage,class TNLabelImage>
EnergyFunction<TImage,TNLabelImage>::EnergyFunction(float alpha, float beta, float gamma):
	m_Alpha(alpha), m_Beta(beta), m_Gamma(gamma)
{
}

template<class TImage,class TNLabelImage>
EnergyFunction<TImage,TNLabelImage>::~EnergyFunction(void)
{
}

template<class TImage,class TNLabelImage>
void
EnergyFunction<TImage,TNLabelImage >::operator()(const ImagePixelType newCenterPixel,
	const NeighborhoodIterator<TImage > &dispIter,
	const NeighborhoodIterator<TNLabelImage > &labelIter,
	const NeighborhoodIterator<TImage > &obsIter,
	const NeighborhoodIterator<TImage > &confidenceIter,
	typename EnergyValueType &combinedEnergy,
	typename EnergyValueType &priorEnergy,
	typename EnergyValueType &observationEnergy) const {

	// save values for debugging
	//m_PriorEnergy = 
	//	this->PriorEnergy(newCenterPixel, mechanicalProperties, labels, it);
	//m_ObservationEnergy = 
	//	this->ObservationEnergy(newCenterPixel, mechanicalProperties, labels, it);
	//m_TotalEnergy = alpha*m_PriorEnergy + (1-alpha)*m_ObservationEnergy;

	// Utotal = alpha*Uprior + (1-alpha)*Uobservation
	priorEnergy = 
		/*m_Alpha**/this->PriorEnergy(newCenterPixel, dispIter, labelIter, obsIter, confidenceIter);
	observationEnergy = 
		/*(1-m_Alpha)**/this->ObservationEnergy(newCenterPixel, dispIter, labelIter, obsIter, 
		confidenceIter);
	combinedEnergy = priorEnergy + observationEnergy;
	//return alpha*this->PriorEnergy(newCenterPixel, mechanicalProperties, labels, it) + 
	//	(1-alpha)*this->ObservationEnergy(newCenterPixel, mechanicalProperties, labels, it);
}

}// end namespace itk
#endif