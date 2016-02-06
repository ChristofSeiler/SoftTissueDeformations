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

#ifndef __itkEnergyFunction_h
#define __itkEnergyFunction_h

#include "itkNeighborhoodIterator.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"
#include "itkImageRegionIterator.h"
#include "MechanicalProperties.h"

using namespace dvfRegularization;
using namespace std;

namespace itk {

template <class TImage,class TNLabelImage>
class ITK_EXPORT EnergyFunction
{
public:
	EnergyFunction();
	EnergyFunction(float alpha, float beta, float gamma);
	~EnergyFunction();

	/** Standard typedefs */
	typedef EnergyFunction Self;
	typedef typename TImage::PixelType ImagePixelType;
	typedef float EnergyValueType;

	/** Capture some typedefs from the template parameters. */
	itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

	/** Reference oeprator. */
	void operator()(const ImagePixelType newCenterPixel,
		const NeighborhoodIterator<TImage > &dispIter,
		const NeighborhoodIterator<TNLabelImage > &labelIter,
		const NeighborhoodIterator<TImage > &obsIter,
		const NeighborhoodIterator<TImage > &confidenceIter,
		typename EnergyValueType &combinedEnergy,
		typename EnergyValueType &priorEnergy,
		typename EnergyValueType &observationEnergy) const;

protected:
	/** Abstract functions */
	virtual EnergyValueType PriorEnergy(const ImagePixelType newCenterPixel,
		const NeighborhoodIterator<TImage > &dispIter,
		const NeighborhoodIterator<TNLabelImage > &labelIter,
		const NeighborhoodIterator<TImage > &obsIter,
		const NeighborhoodIterator<TImage > &confidenceIter) const = 0;
	virtual EnergyValueType ObservationEnergy(const ImagePixelType newCenterPixel,
		const NeighborhoodIterator<TImage > &dispIter,
		const NeighborhoodIterator<TNLabelImage > &labelIter,
		const NeighborhoodIterator<TImage > &obsIter,
		const NeighborhoodIterator<TImage > &confidenceIter) const = 0;
	float m_Alpha;
	float m_Beta;
	float m_Gamma;
};

} // end namespace itk

/*
#define ITK_TEMPLATE_EnergyFunction(_, EXPORT, x, y) \
  namespace itk \
  { \
  _(1(class EXPORT EnergyFunction< ITK_TEMPLATE_1 x >)) \
  namespace Templates { typedef EnergyFunction< ITK_TEMPLATE_1 x > EnergyFunction##y; }\
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkEnergyFunction+-.h"
#endif
*/

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEnergyFunction.txx"
#endif

#endif
