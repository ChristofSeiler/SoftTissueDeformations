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

#include "RegularizationCostFunction.h"
#include "MechanicalProperties.h"
#include <itkEnergyFunction.h>
#include <vtkImageData.h>
#include <limits>
#include "itkMRFRegularizationFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkYoungPoissonEnergyFunction.h"
#include "itkConstantBoundaryCondition.h"

using namespace dvfRegularization;

RegularizationCostFunction::MeasureType RegularizationCostFunction::GetValue( 
	const ParametersType & parameters ) const {

	// run markov own itk filter
	typedef itk::ConstantBoundaryCondition<VectorImageType> BC;
	//typedef itk::ZeroFluxNeumannBoundaryCondition<ImageType> BC;
	typedef itk::MRFRegularizationFilter<VectorImageType,VectorImageType,
		itk::YoungPoissonEnergyFunction<VectorImageType,BC>, BC> MRFFilter;
	MRFFilter::Pointer mrfFilter = MRFFilter::New();
	mrfFilter->SetInput(m_Initial);
	mrfFilter->SetMaximumNumberOfIterations(50);
	mrfFilter->SetDuringIterationUpdate(false);
	mrfFilter->SetStopCriteria(0.01);
	mrfFilter->SetLabels(m_Label);
	mrfFilter->SetMechanicalProperties(m_MechanicalProperties);
	ImageType::Pointer tmpA[3];
	for(unsigned int i = 0; i < 3; i++) tmpA[i] = (*m_ElasticVector)[i];
	itk::YoungPoissonEnergyFunction<VectorImageType,BC> energyFunction(parameters[0], parameters[1], 
		parameters[2], parameters[3], tmpA);
	mrfFilter->SetEnergyFunction(&energyFunction);
	mrfFilter->SetVisualize(false);
	mrfFilter->Update();

	RegularizationCostFunction::MeasureType value 
		= Utils::leastSquare(m_Abaqus, Utils::sampleDVF(mrfFilter->GetOutput()));
	// compare sampled mrf and abaqus images
	cout << "new least square difference = " << value << endl;

	return value;
}

void RegularizationCostFunction::GetDerivative( const ParametersType & parameters,
				   DerivativeType & derivative ) const {
	// no need to implement this at the moment, no used by simplex optimization method
}
