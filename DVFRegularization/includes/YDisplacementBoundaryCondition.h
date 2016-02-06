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
class YDisplacementBoundaryCondition
  : public ImageBoundaryCondition<TImage>
{
public:
  /** Self & superclass typedefs */ 
  typedef YDisplacementBoundaryCondition      Self;
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
  YDisplacementBoundaryCondition()
  { 
	  m_Constant = NumericTraits<PixelType>::Zero;
	  m_MaxDeformationX = NumericTraits<PixelType>::Zero;
	  m_MaxDeformationY = NumericTraits<PixelType>::Zero;
  }

  /** Computes and returns appropriate out-of-bounds values from
   * neighborhood iterator data. */
  virtual PixelType operator()(const OffsetType& point_index,
                               const OffsetType& boundary_offset,
                               const NeighborhoodType *data) const  
  { 
	int linear_index = 0;
	// Return the value of the pixel at the closest boundary point.
	for (unsigned int i = 0; i < ImageDimension; ++i)
    {
		linear_index += (point_index[i] + boundary_offset[i]) * data->GetStride(i);
    }

	  if(boundary_offset[0] < 0) return m_MaxDeformationX;
	  else if(boundary_offset[1] < 0) 
		  return *(reinterpret_cast< PixelType *>( (data->operator[](linear_index)) ));
	  else if(boundary_offset[1] > 0) return -m_MaxDeformationY;
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

  void SetMaxDeformationX(const PixelType &maxDeformationX)
  { m_MaxDeformationX = maxDeformationX; }

  const PixelType &GetMaxDeformationX() const
  { return m_MaxDeformationX; }

  void SetMaxDeformationY(const PixelType &maxDeformationY)
  { m_MaxDeformationY = maxDeformationY; }

  const PixelType &GetMaxDeformationY() const
  { return m_MaxDeformationY; }
  
  /** Tell if the boundary condition can index to any location within
    * the associated iterator's neighborhood or if it has some limited
    * subset (such as none) that it relies upon. */
  bool RequiresCompleteNeighborhood() { return false; }

private:
  PixelType m_Constant;
  PixelType m_MaxDeformationX;
  PixelType m_MaxDeformationY;
};