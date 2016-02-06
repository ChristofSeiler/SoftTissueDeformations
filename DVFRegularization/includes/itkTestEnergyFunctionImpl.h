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

#ifndef __itkTestEnergyFunctionImpl_h
#define __itkTestEnergyFunctionImpl_h

#include "itkEnergyFunction.h"

namespace itk {

template <class TImage, class TBoundaryCondition>
class ITK_EXPORT TestEnergyFunctionImpl : public EnergyFunction<TImage,TBoundaryCondition>
{
public:
	/** Standard class typedefs. */
	typedef TestEnergyFunctionImpl Self;
protected:
	ImagePixelType PriorEnergy(const ImagePixelType newCenterPixel, 
		const vector<MechanicalProperties*>* mechanicalProperties,
		const ImageRegionIterator<LabelType> &labels,
		const NeighborhoodIterator<TImage> &it) const;
	ImagePixelType ObservationEnergy(const ImagePixelType newCenterPixel,
		const vector<MechanicalProperties*>* mechanicalProperties,
		const ImageRegionIterator<LabelType> &labels,
		const NeighborhoodIterator<TImage> &it) const;
};

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTestEnergyFunctionImpl.txx"
#endif

#endif
