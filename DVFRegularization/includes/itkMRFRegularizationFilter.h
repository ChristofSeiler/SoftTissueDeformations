/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Regularize vector displacement field.
 */

#ifndef __itkMRFRegularizationFilter_h
#define __itkMRFRegularizationFilter_h

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkEnergyFunction.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

using namespace dvfRegularization;

namespace itk
{
template <class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
class ITK_EXPORT MRFRegularizationFilter :
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
  typedef float EnergyValueType;

  /** Standard class typedefs. */
  typedef MRFRegularizationFilter Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MRFRegularizationFilter, ImageToImageFilter);
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;

  /** Set/Get the number of iteration of the Iterated Conditional Mode
   * (ICM) algorithm. A default value is set at 50 iterations. */
  itkSetMacro(MaximumNumberOfIterations, unsigned int);

  itkSetMacro(ObservationImage, typename OutputImageType::Pointer);

  itkSetMacro(ConfidenceImage, typename OutputImageType::Pointer);

  itkSetMacro(Labels, typename TNLabelImage::Pointer);

  itkSetMacro(DuringIterationUpdate, bool);

  itkSetMacro(RandomIteration, bool);

  itkSetMacro(StopCriteria, float);

  itkSetMacro(Visualize, bool);
  
  void SetEnergyFunction(typename TEnergyFunction* energyFunction) {
	  m_EnergyFunction = energyFunction;
  }

  /* Get macro for number of iterations */
  itkGetConstReferenceMacro( NumberOfIterations, unsigned int );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimension,
    (itk::Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension), 
		itkGetStaticConstMacro(OutputImageDimension)>));
  /** End concept checking */
#endif

protected:
  MRFRegularizationFilter();
  virtual ~MRFRegularizationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
  virtual void MinimizeFunctional(typename OutputImageType::Pointer output,
	  typename NeighborhoodIterator<OutputImageType>::RadiusType radius);
  virtual float ComputeRMSDFieldEnergy(typename OutputImageType::Pointer output, 
	  typename OutputImageType::Pointer backUp,
	  typename NeighborhoodIterator<OutputImageType>::RadiusType radius);
  virtual float ComputeRMSDDisplacement(typename ImageRegionIterator<OutputImageType> &backUpIter, 
	  typename ImageRegionIterator<OutputImageType> &currentIter);
  virtual void GenerateData();

private:
  MRFRegularizationFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int m_MaximumNumberOfIterations;
  unsigned int m_NumberOfIterations;
  bool m_DuringIterationUpdate;
  float m_StopCriteria;
  bool m_Visualize;
  bool m_RandomIteration;

  typename OutputImageType::Pointer m_Shadow;
  typename OutputImageType::Pointer m_BackupImage;
  typename OutputImageType::Pointer m_ObservationImage;
  typename OutputImageType::Pointer m_ConfidenceImage;
  typename TNLabelImage::Pointer m_Labels;
  TEnergyFunction* m_EnergyFunction;

#ifdef _DEBUG
  ofstream energyFile;
#endif
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRFRegularizationFilter.txx"
#endif

#endif
