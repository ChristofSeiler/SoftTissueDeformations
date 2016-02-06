/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	GPU implemenation of itkMRFRegularizationFilter.
 */

#ifndef __itkMRFRegularizationFilterGPU_h
#define __itkMRFRegularizationFilterGPU_h

#include "itkImageToImageFilter.h"

#include <GL/glew.h>
#include <GL/glut.h>
#include <Cg/cgGL.h>

using namespace dvfRegularization;

namespace itk
{

// Cg vars
CGprofile fragmentProfile;
CGprogram fragmentProgram;
CGparameter image, observation, confidence, dimension;
GLuint glutWindowHandle;
CGcontext cgContext;
// struct for variable parts of GL calls (texture format, float format etc)
struct struct_textureParameters {
	char* name;
	GLenum texTarget;
	GLenum texInternalFormat;
	GLenum texFormat;
	char* shader_source;
}	rect_ati_rgba_32,	// texture rectangles, ATI_texture_float, RGBA, 32 bits
	rect_nv_rgba_32;	// texture rectangles, NV_float_buffer, RGBA, 32 bits

void initGPU();
void cleanUpGPU();

template <class TInputImage, class TOutputImage, class TNLabelImage>
class ITK_EXPORT MRFRegularizationFilterGPU :
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
  typedef TNLabelImage LabelImageType;
  typedef float EnergyValueType;

  /** Standard class typedefs. */
  typedef MRFRegularizationFilterGPU Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MRFRegularizationFilterGPU, ImageToImageFilter);
  
  /** Image typedef support. */
  typedef typename InputImageType::PixelType InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  
  typedef typename InputImageType::RegionType InputImageRegionType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename InputImageType::SizeType InputSizeType;

  /** Set/Get the number of iteration of the Iterated Conditional Mode
   * (ICM) algorithm. A default value is set at 40 iterations. */
  itkSetMacro(NumberOfIterations, unsigned int);

  itkSetMacro(GPUType, unsigned int);

  itkSetMacro(ObservationImage, typename OutputImageType::Pointer);

  itkSetMacro(ConfidenceImage, typename OutputImageType::Pointer);

  itkSetMacro(Labels, typename LabelImageType::Pointer);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimension,
    (itk::Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension), 
		itkGetStaticConstMacro(OutputImageDimension)>));
  /** End concept checking */
#endif

protected:
  MRFRegularizationFilterGPU();
  virtual ~MRFRegularizationFilterGPU();
  void PrintSelf(std::ostream& os, Indent indent) const;
  virtual void GenerateData();

private:
  MRFRegularizationFilterGPU(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int m_NumberOfIterations;
  unsigned int m_GPUType;

  typename OutputImageType::Pointer m_ObservationImage;
  typename OutputImageType::Pointer m_ConfidenceImage;
  typename LabelImageType::Pointer m_Labels;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMRFRegularizationFilterGPU.txx"
#endif

#endif
