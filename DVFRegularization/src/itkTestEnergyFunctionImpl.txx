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

#ifndef _itkTestEnergyFunctionImpl_txx
#define _itkTestEnergyFunctionImpl_txx
#include "itkTestEnergyFunctionImpl.h"

#include "itkEnergyFunction.h"

namespace itk {

template<class TImage, class TBoundaryCondition>
typename TestEnergyFunctionImpl<TImage,TBoundaryCondition>::ImagePixelType 
TestEnergyFunctionImpl<TImage,TBoundaryCondition>::PriorEnergy(
	const ImagePixelType newCenterPixel,
	const vector<MechanicalProperties*>* 
		mechanicalProperties,
	const ImageRegionIterator<LabelType> &labels,
	const NeighborhoodIterator<TImage> &it) const
{
	ImagePixelType priorEnergy = NumericTraits<ImagePixelType>::Zero;

	// set weight with the help of mechanical properties
	float youngsModulus = (*mechanicalProperties)[labels.Get()]->getYoungsModulus();
	float poissonRatio = (*mechanicalProperties)[labels.Get()]->getPoissonRatio();

	for(unsigned int i = 0; i < 9; i++) {
		ImagePixelType distance = (newCenterPixel - it.GetPixel(i));
		priorEnergy += (distance * distance);
	}

	return priorEnergy;
}

template<class TImage, class TBoundaryCondition>
typename TestEnergyFunctionImpl<TImage,TBoundaryCondition>::ImagePixelType 
TestEnergyFunctionImpl<TImage,TBoundaryCondition>::ObservationEnergy(
	const ImagePixelType newCenterPixel,
	const vector<MechanicalProperties*>* 
		mechanicalProperties,
	const ImageRegionIterator<LabelType> &labels,
	const NeighborhoodIterator<TImage> &it) const
{
	ImagePixelType observationEnergy = NumericTraits<ImagePixelType>::Zero;
	return observationEnergy;
}

}// end namespace itk
#endif
