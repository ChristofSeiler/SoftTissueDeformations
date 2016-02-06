/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	
 */

#ifndef __itkInterpolateNextStepFilter_h
#define __itkInterpolateNextStepFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage, class TInterpolator>
class ITK_EXPORT InterpolateNextStepFilter :
    public ImageToImageFilter< TInputImage, TOutputImage>
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
  typedef itk::Image<float,2> FloatImageType;
  typedef TInterpolator InterpolatorType;

  /** Standard class typedefs. */
  typedef InterpolateNextStepFilter Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(InterpolateNextStepFilter, ImageToImageFilter);
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  typedef typename InputImageType::SizeType InputSizeType;

  itkSetMacro(NumberOfPixels, unsigned int);

protected:
  InterpolateNextStepFilter();
  virtual ~InterpolateNextStepFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
  virtual void GenerateData();

private:
  InterpolateNextStepFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int m_NumberOfPixels;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInterpolateNextStepFilter.txx"
#endif

#endif
