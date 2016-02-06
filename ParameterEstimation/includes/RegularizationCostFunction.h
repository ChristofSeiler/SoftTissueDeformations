/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Implementation of Metropolis Algorithm for sampling the associated Gibbs 
 *			distribution and minimizing the energy function.
 */

#pragma once

#include "itkSingleValuedCostFunction.h"
#include "Definitions.h"

namespace dvfRegularization {
	class MechanicalProperties;
  
/** \class SingleValuedCostFunction
 * \brief This class is a base for the CostFunctions returning a 
 * single value
 *
 * \ingroup Numerics Optimizers
 */
class RegularizationCostFunction : public itk::SingleValuedCostFunction 
{
public:
  typedef RegularizationCostFunction   Self;
  typedef itk::SingleValuedCostFunction     Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;
  itkNewMacro( Self );
  itkTypeMacro( RegularizationCostFunction, SingleValuedCostFunction );

  enum { SpaceDimension=4 };

  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef Superclass::MeasureType                 MeasureType;

  void SetInitial(VectorImageType* initial) { m_Initial = initial; }
  void SetElasticVector(std::vector<ImageType::Pointer>* elasticVector) { m_ElasticVector = elasticVector; }
  void SetLabel(ImageType* label) { m_Label = label; }
  void SetMechanicalProperties(std::vector<MechanicalProperties*>* mechanicalProperties) {
	  m_MechanicalProperties = mechanicalProperties; }
  void SetAbaqus(VectorImageType* abaqus) { m_Abaqus = abaqus; }

  /** This method returns the value of the cost function corresponding
    * to the specified parameters.    */ 
  MeasureType GetValue( const ParametersType & parameters ) const;

  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters.   */ 
  void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const;
  
  unsigned int GetNumberOfParameters(void) const {
    return SpaceDimension;
  }

protected:
  RegularizationCostFunction() {};
  virtual ~RegularizationCostFunction() {};

private:
  RegularizationCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  VectorImageType* m_Initial;
  std::vector<ImageType::Pointer>* m_ElasticVector;
  ImageType* m_Label;
  std::vector<MechanicalProperties*>* m_MechanicalProperties;
  VectorImageType* m_Abaqus;

};

} // end namespace dvfRegularization
