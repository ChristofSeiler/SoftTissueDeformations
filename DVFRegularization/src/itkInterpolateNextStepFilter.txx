/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	
 */

#ifndef _itkInterpolateNextStepFilter_txx
#define _itkInterpolateNextStepFilter_txx
#include "itkInterpolateNextStepFilter.h"

#include "itkResampleImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkConstantBoundaryCondition.h"

namespace itk
{

template<class TInputImage, class TOutputImage, class TInterpolator>
InterpolateNextStepFilter<TInputImage,TOutputImage,TInterpolator>
::InterpolateNextStepFilter(void):
  m_NumberOfPixels(0)
{
}

template<class TInputImage, class TOutputImage, class TInterpolator>
InterpolateNextStepFilter<TInputImage,TOutputImage,TInterpolator>
::~InterpolateNextStepFilter(void)
{
}

template<class TInputImage, class TOutputImage, class TInterpolator>
void
InterpolateNextStepFilter<TInputImage,TOutputImage,TInterpolator>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent <<" Resampling to next step " << std::endl;

  os << indent <<" Number of pixels: " << 
    m_NumberOfPixels << std::endl;

}// end PrintSelf

template<class TInputImage, class TOutputImage, class TInterpolator>
void
InterpolateNextStepFilter<TInputImage,TOutputImage,TInterpolator>
::GenerateData()
{
	// set input image to output image
	InputImageType::ConstPointer input = this->GetInput();

	// split
	FloatImageType::Pointer single[2];
	FloatImageType::Pointer singleOut[2];
	single[0] = FloatImageType::New();
	single[1] = FloatImageType::New();
	FloatImageType::IndexType start;
	FloatImageType::SizeType size;
	size[0] = input->GetRequestedRegion().GetSize()[0];
	size[1] = input->GetRequestedRegion().GetSize()[1];
	start[0] = 0; start[1] = 0;
	FloatImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	single[0]->SetRegions(region);
	single[1]->SetRegions(region);
	single[0]->Allocate();
	single[1]->Allocate();
	itk::ImageRegionIterator<FloatImageType> floatIter1(single[0], single[0]->GetRequestedRegion());
	itk::ImageRegionIterator<FloatImageType> floatIter2(single[1], single[1]->GetRequestedRegion());
	itk::ImageRegionConstIterator<InputImageType> imageIter(input, input->GetRequestedRegion());
	for(imageIter.GoToBegin(), floatIter1.GoToBegin(), floatIter2.GoToBegin(); 
		!imageIter.IsAtEnd(); ++imageIter, ++floatIter1, ++floatIter2) {
			floatIter1.Set(imageIter.Get()[0]);
			floatIter2.Set(imageIter.Get()[1]);
	}
	// resample
	for(int i = 0; i < 2; ++i) {
		typedef itk::ResampleImageFilter<FloatImageType,FloatImageType> ResampleFilterType;
		ResampleFilterType::Pointer filter = ResampleFilterType::New();
		typedef itk::AffineTransform< double, 2 > TransformType;
		TransformType::Pointer transform = TransformType::New();
		transform->Scale((input->GetBufferedRegion().GetSize()[0]-1.0)/(std::sqrt((float)m_NumberOfPixels)-1.0));
		filter->SetTransform( transform );
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		filter->SetInterpolator( interpolator );
		double spacing[ 2 ];
		spacing[0] = 1.0; spacing[1] = 1.0;
		filter->SetOutputSpacing( spacing );
		double origin[ 2 ];
		origin[0] = 0.0; origin[1] = 0.0;
		filter->SetOutputOrigin( origin );
		FloatImageType::SizeType size;
		size[0] = std::sqrt((float)m_NumberOfPixels); size[1] = std::sqrt((float)m_NumberOfPixels);
		filter->SetSize( size );
		filter->SetInput( single[i] );
		filter->Update();
		singleOut[i] = filter->GetOutput();
	}
	// combine
	OutputImageType::Pointer image = this->GetOutput();
	OutputImageType::IndexType newStart;
	OutputImageType::SizeType newSize;
	newSize.Fill(1);
	newSize[0] = singleOut[0]->GetRequestedRegion().GetSize()[0];
	newSize[1] = singleOut[0]->GetRequestedRegion().GetSize()[1];
	newStart.Fill(0);
	OutputImageType::RegionType newRegion;
	newRegion.SetSize( newSize );
	newRegion.SetIndex( newStart );
	image->SetRegions( newRegion );
	image->Allocate();
	itk::ImageRegionIterator<FloatImageType> floatIter1After(singleOut[0], singleOut[0]->GetRequestedRegion());
	itk::ImageRegionIterator<FloatImageType> floatIter2After(singleOut[1], singleOut[1]->GetRequestedRegion());
	itk::ImageRegionIterator<OutputImageType> imageIterNew(image, image->GetRequestedRegion());
	for(imageIter.GoToBegin(), floatIter1After.GoToBegin(), floatIter2After.GoToBegin(); 
		!imageIterNew.IsAtEnd(); ++imageIterNew, ++floatIter1After, ++floatIter2After) {
			OutputImageType::PixelType value;
			value.Fill(0.0);
			value[0] = floatIter1After.Get();
			value[1] = floatIter2After.Get();
			imageIterNew.Set(value);
	}

}// end GenerateData

} // end namespace itk

#endif