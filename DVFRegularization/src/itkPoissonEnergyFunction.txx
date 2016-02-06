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

#ifndef _itkPoissonEnergyFunction_txx
#define _itkPoissonEnergyFunction_txx
#include "itkPoissonEnergyFunction.h"

#include "itkEnergyFunction.h"
#include "itkMedianImageFilter.h"

namespace itk {

template<class TImage,class TBoundaryCondition>
PoissonEnergyFunction<TImage,TBoundaryCondition>::PoissonEnergyFunction(
	float alpha, float beta0, float beta1, float gamma, MedianImageType::Pointer medianImages[3]) : 
	EnergyFunction<TImage,TBoundaryCondition>(alpha),
	m_Beta0(beta0),
	m_Beta1(beta1),
	m_Gamma(gamma)
{
}

template<class TImage,class TBoundaryCondition>
PoissonEnergyFunction<TImage,TBoundaryCondition>::~PoissonEnergyFunction() 
{
}

template<class TImage,class TBoundaryCondition>
typename PoissonEnergyFunction<TImage,TBoundaryCondition>::EnergyValueType 
PoissonEnergyFunction<TImage,TBoundaryCondition>::PriorEnergy(
	const ImagePixelType newCenterPixel,
	const vector<MechanicalProperties*>* mechanicalProperties,
	const NeighborhoodIterator<TImage,TBoundaryCondition> &dispIter,
	const NeighborhoodIterator<TImage,ZeroFluxNeumannBoundaryCondition<TImage> > &labelIter) const 
{
	EnergyValueType priorEnergy = 0;

	// get displacement in x
	// ratio
	float ratio = labelIter.GetCenterPixel()[1]/labelIter.GetPixel(7)[1];

	// first derivatives
	float delta1 = newCenterPixel[1] - dispIter.GetPixel(1)[1];
	float delta2 = dispIter.GetPixel(7)[1] - newCenterPixel[1];
	
	// compare
	float diff = abs(ratio * delta1 - delta2);

	// add penalty for high displacements

	// update energy
	priorEnergy += diff;

	return priorEnergy;
}

template<class TImage,class TBoundaryCondition>
typename PoissonEnergyFunction<TImage,TBoundaryCondition>::EnergyValueType 
PoissonEnergyFunction<TImage,TBoundaryCondition>::ObservationEnergy(
	const ImagePixelType newCenterPixel,
	const vector<MechanicalProperties*>* mechanicalProperties,
	const NeighborhoodIterator<TImage,TBoundaryCondition> &dispIter,
	const NeighborhoodIterator<TImage,ZeroFluxNeumannBoundaryCondition<TImage> > &labelIter) const 
{
	EnergyValueType observationEnergy = 0;
	return observationEnergy;
}

}// end namespace itk
#endif
