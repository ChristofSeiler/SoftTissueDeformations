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

#ifndef _itkYoungPoissonEnergyFunction_txx
#define _itkYoungPoissonEnergyFunction_txx
#include "itkYoungPoissonEnergyFunction.h"

#include "itkEnergyFunction.h"
#include "itkMedianImageFilter.h"

namespace itk {

template<class TImage,class TNLabelImage>
YoungPoissonEnergyFunction<TImage,TNLabelImage>::YoungPoissonEnergyFunction(
	float alpha, float beta, float gamma) : 
	EnergyFunction<TImage,TNLabelImage>(alpha,beta,gamma)
{
}

template<class TImage,class TNLabelImage>
YoungPoissonEnergyFunction<TImage,TNLabelImage>::~YoungPoissonEnergyFunction() 
{
}

template<class TImage,class TNLabelImage>
typename YoungPoissonEnergyFunction<TImage,TNLabelImage>::EnergyValueType 
YoungPoissonEnergyFunction<TImage,TNLabelImage>::PriorEnergy(
	const ImagePixelType newCenterPixel,
	const NeighborhoodIterator<TImage > &dispIter,
	const NeighborhoodIterator<TNLabelImage > &labelIter,
	const NeighborhoodIterator<TImage > &obsIter,
	const NeighborhoodIterator<TImage > &confidenceIter) const 
{
	EnergyValueType priorEnergy = 0;

	// get ym information of center pixel
	float lCenter = labelIter.GetPixel(4);
	
	// ratio in x left
	float ratio = lCenter/labelIter.GetPixel(3);
	// first derivatives in x
	float delta1 = newCenterPixel[0] - dispIter.GetPixel(5)[0];
	float delta2 = dispIter.GetPixel(3)[0] - newCenterPixel[0];
	// compare
	float diff = ratio * delta1 - delta2;
	// update energy
	priorEnergy += diff*diff;

	// ratio in x right
	ratio = lCenter/labelIter.GetPixel(5);
	// first derivatives in x
	delta1 = newCenterPixel[0] - dispIter.GetPixel(3)[0];
	delta2 = dispIter.GetPixel(5)[0] - newCenterPixel[0];
	// compare
	diff = ratio * delta1 - delta2;
	// update energy
	priorEnergy += diff*diff;

	// ratio in y top
	ratio = lCenter/labelIter.GetPixel(1);
	// first derivatives in y
	delta1 = newCenterPixel[1] - dispIter.GetPixel(7)[1];
	delta2 = dispIter.GetPixel(1)[1] - newCenterPixel[1];
	// compare
	diff = ratio * delta1 - delta2;
	// update energy
	priorEnergy += diff*diff;

	// ratio in y bottom
	ratio = lCenter/labelIter.GetPixel(7);
	// first derivatives in y
	delta1 = newCenterPixel[1] - dispIter.GetPixel(1)[1];
	delta2 = dispIter.GetPixel(7)[1] - newCenterPixel[1];
	// compare
	diff = ratio * delta1 - delta2;
	// update energy
	priorEnergy += diff*diff;

	// ratio in north east
	float deltav1[2];
	float deltav2[2];
	float diffv[2];
	ratio = lCenter/labelIter.GetPixel(2);
	deltav1[0] = newCenterPixel[0] - dispIter.GetPixel(6)[0];
	deltav1[1] = newCenterPixel[1] - dispIter.GetPixel(6)[1];
	deltav2[0] = dispIter.GetPixel(2)[0] - newCenterPixel[0];
	deltav2[1] = dispIter.GetPixel(2)[1] - newCenterPixel[1];
	diffv[0] = ratio * deltav1[0] - deltav2[0];
	diffv[1] = ratio * deltav1[1] - deltav2[1];
	priorEnergy += diffv[0]*diffv[0]+diffv[1]*diffv[1];

	// ratio in north west
	ratio = lCenter/labelIter.GetPixel(0);
	deltav1[0] = newCenterPixel[0] - dispIter.GetPixel(8)[0];
	deltav1[1] = newCenterPixel[1] - dispIter.GetPixel(8)[1];
	deltav2[0] = dispIter.GetPixel(0)[0] - newCenterPixel[0];
	deltav2[1] = dispIter.GetPixel(0)[1] - newCenterPixel[1];
	diffv[0] = ratio * deltav1[0] - deltav2[0];
	diffv[1] = ratio * deltav1[1] - deltav2[1];
	priorEnergy += diffv[0]*diffv[0]+diffv[1]*diffv[1];
	
	// ratio in south west
	ratio = lCenter/labelIter.GetPixel(6);
	deltav1[0] = newCenterPixel[0] - dispIter.GetPixel(2)[0];
	deltav1[1] = newCenterPixel[1] - dispIter.GetPixel(2)[1];
	deltav2[0] = dispIter.GetPixel(6)[0] - newCenterPixel[0];
	deltav2[1] = dispIter.GetPixel(6)[1] - newCenterPixel[1];
	diffv[0] = ratio * deltav1[0] - deltav2[0];
	diffv[1] = ratio * deltav1[1] - deltav2[1];
	priorEnergy += diffv[0]*diffv[0]+diffv[1]*diffv[1];
	
	// ratio in south east
	ratio = lCenter/labelIter.GetPixel(8);
	deltav1[0] = newCenterPixel[0] - dispIter.GetPixel(0)[0];
	deltav1[1] = newCenterPixel[1] - dispIter.GetPixel(0)[1];
	deltav2[0] = dispIter.GetPixel(8)[0] - newCenterPixel[0];
	deltav2[1] = dispIter.GetPixel(8)[1] - newCenterPixel[1];
	diffv[0] = ratio * deltav1[0] - deltav2[0];
	diffv[1] = ratio * deltav1[1] - deltav2[1];
	priorEnergy += diffv[0]*diffv[0]+diffv[1]*diffv[1];

	return priorEnergy;
}

template<class TImage,class TNLabelImage>
typename YoungPoissonEnergyFunction<TImage,TNLabelImage>::EnergyValueType 
YoungPoissonEnergyFunction<TImage,TNLabelImage>::ObservationEnergy(
	const ImagePixelType newCenterPixel,
	const NeighborhoodIterator<TImage > &dispIter,
	const NeighborhoodIterator<TNLabelImage > &labelIter,
	const NeighborhoodIterator<TImage > &obsIter,
	const NeighborhoodIterator<TImage > &confidenceIter) const 
{
	EnergyValueType observationEnergy = 0;

	// PCA or Median Root Prior?

	/*
	for(unsigned int j = 0; j < 2; j++) {
		EnergyValueType diff = abs(m_MedianImages[j]->GetPixel(it.GetIndex()) - newCenterPixel[j]);
		observationEnergy += m_Gamma * diff;
	}
	*/

	/*
	// Mean Root Prior used as Mean Root Observation
	// The penalty of MRP is set according to how much the center pixel differs
	// from the local median.
	// alpha between 0.345106288 and 0.345106289

	// monotonic incresing/decreasing
	for(unsigned int j = 0; j < 2; j++) {
		// We have to copy the pixels so we can run std::nth_element.
		std::vector<EnergyValueType> pixels;
		typename std::vector<EnergyValueType>::iterator medianIterator;

		// Walk the neighborhood
		for (unsigned int i = 0; i < dispIter.Size(); ++i)
		{
			pixels.push_back( dispIter.GetPixel(i)[j] );
		}

		// Get the median value
		unsigned int medianPosition = dispIter.Size() / 2;
		medianIterator = pixels.begin() + medianPosition;
		std::nth_element(pixels.begin(), medianIterator, pixels.end());

		EnergyValueType diff = abs(*medianIterator - newCenterPixel[j]);
		//EnergyValueType normalized = diff/abs(newCenterPixel[j]);
		//EnergyValueType gamma = 1*abs(newCenterPixel[j]);
		observationEnergy += m_Beta0 * diff;
	}*/

	/*
	// mean
	TImage::PixelType mean;
	mean[0] = 0; mean[1] = 0; mean[2] = 0;
	unsigned int size = dispIter.Size();
	for (unsigned int i = 0; i < size; ++i) {
		mean += dispIter.GetPixel(i);
	}
	mean *= 1.0/size;
	for(unsigned int j = 0; j < 2; j++)
		observationEnergy += abs(mean[j] - newCenterPixel[j]);
	*/

	if(confidenceIter.GetCenterPixel()[0] == 1)
		observationEnergy += 1000*abs(newCenterPixel[0] - obsIter.GetCenterPixel()[0]);
	if(confidenceIter.GetCenterPixel()[1] == 1)
		observationEnergy += 1000*abs(newCenterPixel[1] - obsIter.GetCenterPixel()[1]);

	return observationEnergy;
}

}// end namespace itk
#endif
