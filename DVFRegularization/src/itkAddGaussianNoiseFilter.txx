/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Add gaussian noise to an image.
 */

#ifndef _itkAddGaussianNoiseFilter_txx
#define _itkAddGaussianNoiseFilter_txx
#include "itkAddGaussianNoiseFilter.h"

#include "itkNormalVariateGenerator.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
AddGaussianNoiseFilter<TInputImage, TOutputImage>
::AddGaussianNoiseFilter() : m_Mean(0), m_SD(1)
{
}

template< class TInputImage, class TOutputImage>
void
AddGaussianNoiseFilter< TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
{
	// Allocate output
	typename OutputImageType::Pointer output = this->GetOutput();
	typename  InputImageType::ConstPointer input  = this->GetInput();

	// support progress methods/callbacks
	ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  	// gauss normal distribution test
	typedef itk::Statistics::NormalVariateGenerator GeneratorType;
	GeneratorType::Pointer generator = GeneratorType::New();
	generator->Initialize( (int) 2003 );

	typedef itk::ImageRegionConstIterator<InputImageType> ItTypeInput;
	typedef itk::ImageRegionIterator<OutputImageType> ItTypeOuput;
	ItTypeInput itInput( input, input->GetRequestedRegion() );
	ItTypeOuput itOutput( output, input->GetRequestedRegion() );
	for (itInput = itInput.Begin(); !itInput.IsAtEnd(); ++itInput, ++itOutput) {
		itOutput.Set( itInput.Get() + m_Mean + m_SD * generator->GetVariate() );
		progress.CompletedPixel();
	}
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput>
void
AddGaussianNoiseFilter<TInputImage, TOutput>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Mean: " << m_Mean << std::endl;
  os << indent << "SD: " << m_SD << std::endl;

}

} // end namespace itk

#endif
