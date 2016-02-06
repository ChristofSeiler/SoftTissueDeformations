/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Add gaussian noise to an image.
 */

#ifndef __itkAddGaussianNoiseFilter_h
#define __itkAddGaussianNoiseFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"

namespace itk
{
template <class TInputImage, class TOutputImage >
class ITK_EXPORT AddGaussianNoiseFilter :
    public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage InputImageType;
  typedef TOutputImage OutputImageType;

  /** Standard class typedefs. */
  typedef AddGaussianNoiseFilter Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AddGaussianNoiseFilter, ImageToImageFilter);
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;

  itkSetMacro(Mean, int);
  itkSetMacro(SD, int);

  itkGetConstReferenceMacro(Mean, int);
  itkGetConstReferenceMacro(SD, int);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimension,
    (itk::Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension), 
		itkGetStaticConstMacro(OutputImageDimension)>));
  /** End concept checking */
#endif

protected:
  AddGaussianNoiseFilter();
  virtual ~AddGaussianNoiseFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            int threadId );

private:
  AddGaussianNoiseFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  int m_Mean;
  int m_SD;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAddGaussianNoiseFilter.txx"
#endif

#endif
