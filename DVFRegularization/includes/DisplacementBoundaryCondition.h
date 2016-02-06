/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	
 */

#pragma once

#include "itkNeighborhood.h"
#include "itkNumericTraits.h"
#include "itkImageBoundaryCondition.h"

/** \class ConstantBoundaryCondition
 * \brief This boundary condition returns a constant value for out-of-bounds
 * image pixels.
 * 
 * For example, invoking this function object with a constant value of zero
 * (the default) on each out-of-bounds element of a 7x5 iterator that masks a
 * region at an image corner 
 * (iterator is centered on the 2): 
 *
 *               * * * * * * * 
 *               * * * * * * *
 *               * * 1 2 3 4 5  (where * denotes pixels that lie 
 *               * * 3 3 5 5 6          outside of the image boundary)
 *               * * 4 4 6 7 8
 *
 * would produce the following neighborhood of values:
 *
 *               0.9 1 0 0 0 0 0
 *               0.9 1 0 0 0 0 0
 *               0.9 1 1 2 3 4 5
 *               0.9 1 3 3 5 5 6
 *               0.9 1 4 4 6 7 8
 *
 * 
 * \note If you are using an image with Array as the pixel type, you will need 
 * to set the constant explicitly with an array of the appropriate length. This 
 * is also true if your image type is a VectorImage.
 * 
 * \sa ImageBoundaryCondition
 *
 * \ingroup DataRepresentation
 * \ingroup ImageObjects
 */

using namespace itk;

template<class TImage>
class DisplacementBoundaryCondition
  : public ImageBoundaryCondition<TImage>
{
public:
  /** Self & superclass typedefs */ 
  typedef DisplacementBoundaryCondition      Self;
  typedef ImageBoundaryCondition<TImage> Superclass;
  
  /** Extract information from the image type */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::PixelPointerType PixelPointerType;
  typedef typename Superclass::IndexType        IndexType;
  typedef typename Superclass::OffsetType       OffsetType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
    
  typedef typename Superclass::NeighborhoodAccessorFunctorType 
                                 NeighborhoodAccessorFunctorType;

  /** Save the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Default constructor. */
  DisplacementBoundaryCondition()
  { 
	  m_Constant = NumericTraits<PixelType>::Zero;
	  m_MaxDisplacementRight = NumericTraits<PixelType>::Zero;
	  m_MaxDisplacementLeft = NumericTraits<PixelType>::Zero;
	  m_MaxDisplacementTop = NumericTraits<PixelType>::Zero;
	  m_MaxDisplacementBottom = NumericTraits<PixelType>::Zero;
  }

  /** Computes and returns appropriate out-of-bounds values from
   * neighborhood iterator data. */
  virtual PixelType operator()(const OffsetType& point_index,
                               const OffsetType& boundary_offset,
                               const NeighborhoodType *data) const  
  { 
	  if(boundary_offset[0] < 0) return m_MaxDisplacementRight;
	  else if(boundary_offset[0] > 0) return m_MaxDisplacementLeft;
	  else if(boundary_offset[1] > 0) return m_MaxDisplacementTop;
	  else if(boundary_offset[1] < 0) return m_MaxDisplacementBottom;
	  return m_Constant; 
  }

  /** Computes and returns the appropriate pixel value from
   * neighborhood iterator data, using the functor. */
  virtual PixelType operator()(
      const OffsetType& point_index,
      const OffsetType& boundary_offset,
      const NeighborhoodType *data,
      const NeighborhoodAccessorFunctorType &neighborhoodAccessorFunctor) const  
  { 
	  return this->operator()(point_index, boundary_offset, data); 
  }
  
  /** Set the value of the constant. */
  void SetConstant(const PixelType &c)
    {  m_Constant = c; }

  /** Get the value of the constant. */
  const PixelType &GetConstant() const
    {  return m_Constant;  }

  void SetMaxDisplacementRight(const PixelType &maxDisplacementRight)
  { m_MaxDisplacementRight = maxDisplacementRight; }

  const PixelType &GetMaxDisplacementRight() const
  { return m_MaxDisplacementRight; }

  void SetMaxDisplacementLeft(const PixelType &maxDisplacementLeft)
  { m_MaxDisplacementLeft = maxDisplacementLeft; }

  const PixelType &GetMaxDisplacementLeft() const
  { return m_MaxDisplacementLeft; }

  void SetMaxDisplacementTop(const PixelType &maxDisplacementTop)
  { m_MaxDisplacementTop = maxDisplacementTop; }

  const PixelType &GetMaxDisplacementTop() const
  { return m_MaxDisplacementTop; }

  void SetMaxDisplacementBottom(const PixelType &maxDisplacementBottom)
  { m_MaxDisplacementBottom = maxDisplacementBottom; }

  const PixelType &GetMaxDisplacementBottom() const
  { return m_MaxDisplacementBottom; }
  
  /** Tell if the boundary condition can index to any location within
    * the associated iterator's neighborhood or if it has some limited
    * subset (such as none) that it relies upon. */
  bool RequiresCompleteNeighborhood() { return false; }

private:
  PixelType m_Constant;
  PixelType m_MaxDisplacementTop;
  PixelType m_MaxDisplacementBottom;
  PixelType m_MaxDisplacementLeft;
  PixelType m_MaxDisplacementRight;
};