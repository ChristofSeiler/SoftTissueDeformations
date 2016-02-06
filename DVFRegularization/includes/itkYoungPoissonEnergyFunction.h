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

#ifndef __itkYoungPoissonEnergyFunction_h
#define __itkYoungPoissonEnergyFunction_h

#include "itkEnergyFunction.h"

namespace itk {

template <class TImage,class TNLabelImage>
class ITK_EXPORT YoungPoissonEnergyFunction : 
	public EnergyFunction<TImage,TNLabelImage>
{
public:
	/** Standard class typedefs. */
	typedef YoungPoissonEnergyFunction Self;

	YoungPoissonEnergyFunction(float alpha, float beta, float gamma);
	~YoungPoissonEnergyFunction();

protected:
	EnergyValueType PriorEnergy(const ImagePixelType newCenterPixel, 
		const NeighborhoodIterator<TImage > &dispIter,
		const NeighborhoodIterator<TNLabelImage > &labelIter,
		const NeighborhoodIterator<TImage > &obsIter,
		const NeighborhoodIterator<TImage > &confidenceIter) const;
	EnergyValueType ObservationEnergy(const ImagePixelType newCenterPixel,
		const NeighborhoodIterator<TImage > &dispIter,
		const NeighborhoodIterator<TNLabelImage > &labelIter,
		const NeighborhoodIterator<TImage > &obsIter,
		const NeighborhoodIterator<TImage > &confidenceIter) const;
};

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkYoungPoissonEnergyFunction.txx"
#endif

#endif
