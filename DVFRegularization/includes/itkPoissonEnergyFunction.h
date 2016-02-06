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

#ifndef __itkPoissonEnergyFunction_h
#define __itkPoissonEnergyFunction_h

#include "itkEnergyFunction.h"

namespace itk {

template <class TImage,class TBoundaryCondition>
class ITK_EXPORT PoissonEnergyFunction : 
	public EnergyFunction<TImage,TBoundaryCondition>
{
public:
	/** Standard class typedefs. */
	typedef PoissonEnergyFunction Self;

	typedef itk::Image<float,3> MedianImageType;

	PoissonEnergyFunction(float alpha, float beta0, float beta1,
		float gamma, MedianImageType::Pointer medianImages[3]);
	~PoissonEnergyFunction();

protected:
	EnergyValueType PriorEnergy(const ImagePixelType newCenterPixel, 
		const vector<MechanicalProperties*>* mechanicalProperties,
		const NeighborhoodIterator<TImage,TBoundaryCondition> &dispIter,
		const NeighborhoodIterator<TImage,ZeroFluxNeumannBoundaryCondition<TImage> > &labelIter) const;
	EnergyValueType ObservationEnergy(const ImagePixelType newCenterPixel,
		const vector<MechanicalProperties*>* mechanicalProperties,
		const NeighborhoodIterator<TImage,TBoundaryCondition> &dispIter,
		const NeighborhoodIterator<TImage,ZeroFluxNeumannBoundaryCondition<TImage> > &labelIter) const;
private:
	float m_Beta0;
	float m_Beta1;
	float m_Gamma;
	MedianImageType::Pointer m_MedianImages[3];
};

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPoissonEnergyFunction.txx"
#endif

#endif
