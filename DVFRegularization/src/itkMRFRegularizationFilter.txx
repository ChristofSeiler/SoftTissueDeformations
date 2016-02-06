/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Regularize vector displacement field.
 */

#ifndef _itkMRFRegularizationFilter_txx
#define _itkMRFRegularizationFilter_txx
#include "itkMRFRegularizationFilter.h"

#include "Utils.h"

#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkMinimumDecisionRule.h"

#include <iostream>
#include <fstream>

namespace itk
{
template<class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
MRFRegularizationFilter<TInputImage,TOutputImage,TNLabelImage,TEnergyFunction>
::MRFRegularizationFilter(void):
  m_MaximumNumberOfIterations(50),
  m_NumberOfIterations(0),
  m_DuringIterationUpdate(false),
  m_StopCriteria(0.0),
  m_Visualize(false),
  m_RandomIteration(false)
{

  if( (int)InputImageDimension != (int)OutputImageDimension )
    {
    OStringStream msg;
    msg << "Input image dimension: " << InputImageDimension << 
		" != output image dimension: " 
		<< OutputImageDimension; 
    throw ExceptionObject(__FILE__, __LINE__,msg.str().c_str(),ITK_LOCATION);
    }
}

template<class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
MRFRegularizationFilter<TInputImage,TOutputImage,TNLabelImage,TEnergyFunction>
::~MRFRegularizationFilter(void)
{

}

template<class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
void
MRFRegularizationFilter<TInputImage,TOutputImage,TNLabelImage,TEnergyFunction>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent <<" MRF Regulariziation filter object " << std::endl;

  os << indent <<" Maximum number of iterations: " << 
    m_MaximumNumberOfIterations << std::endl;

  os << indent <<" Number of iterations: " << 
    m_NumberOfIterations << std::endl;

}// end PrintSelf

template<class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
void
MRFRegularizationFilter<TInputImage,TOutputImage,TNLabelImage,TEnergyFunction>
::GenerateData()
{
	// set input image to output image
	InputImageType::ConstPointer input = this->GetInput();
	OutputImageType::Pointer output = this->GetOutput();
	output->SetRegions(input->GetRequestedRegion());
	output->Allocate();
	ImageRegionConstIterator<InputImageType> in(input, input->GetRequestedRegion());
	ImageRegionIterator<OutputImageType> out(output, output->GetRequestedRegion());
	in.GoToBegin(); out.GoToBegin();
	for(in.GoToBegin(), out.GoToBegin();!in.IsAtEnd(); ++out, ++in)
		out.Set(in.Get());

	// init shadow image
	m_Shadow = OutputImageType::New();
	m_Shadow->SetRegions(output->GetRequestedRegion());
	m_Shadow->Allocate();

	// init backup image
	m_BackupImage = OutputImageType::New();
	m_BackupImage->SetRegions(output->GetRequestedRegion());
	m_BackupImage->Allocate();
	ImageRegionIterator<OutputImageType> backupIter(m_BackupImage, m_BackupImage->GetRequestedRegion());

	if(m_Visualize) {
		// todo
	}

	// radius
	NeighborhoodIterator<OutputImageType>::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;

	// iterate to convergence
	//EnergyValueType previousEnergy = ComputeFieldEnergy(output, radius);
	//EnergyValueType currentEnergy = 0;
	EnergyValueType rmsd = 0;
#ifdef _DEBUG
	ofstream energyFile;
	int numberOfPixels = output->GetRequestedRegion().GetSize()[0];
	numberOfPixels *= numberOfPixels;
	std::stringstream stream;
	stream << "energyFile" << numberOfPixels << ".txt";
	energyFile.open (stream.str().c_str());
#endif
	for(; m_NumberOfIterations < m_MaximumNumberOfIterations; ++m_NumberOfIterations) {

		// backup current image
		for(backupIter.GoToBegin(), out.GoToBegin(); !backupIter.IsAtEnd(); ++backupIter, ++out)
			backupIter.Set(out.Get());
		
		// run ICM
		MinimizeFunctional(output, radius);

		if(m_StopCriteria != 0) {
			// did it converge?
			rmsd = ComputeRMSDFieldEnergy(m_BackupImage, output, radius);
			if(rmsd < m_StopCriteria) {
				break;
			}
		}

		// render
		if(m_Visualize) {
			// TODO
		}
		
		// log
#ifdef _DEBUG
		cout << "Finished with iteration " << m_NumberOfIterations << ", current energy = " 
			<< rmsd << endl;
		energyFile << "Finished with iteration " << m_NumberOfIterations << ", current energy = " 
			<< rmsd << endl;
#endif
		//previousEnergy = currentEnergy;
	}
#ifdef _DEBUG
	energyFile.close();
#endif

}// end GenerateData

template<class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
void
MRFRegularizationFilter<TInputImage,TOutputImage,TNLabelImage,TEnergyFunction>
::MinimizeFunctional(typename OutputImageType::Pointer output,
					 typename NeighborhoodIterator<OutputImageType>::RadiusType radius)
{
	// build up neighborhood iteration stuff
	NeighborhoodIterator<OutputImageType> 
		dispIter(radius, output, output->GetRequestedRegion());

	NeighborhoodIterator<TNLabelImage > 
		labelIter(radius, m_Labels, m_Labels->GetRequestedRegion());

	NeighborhoodIterator<OutputImageType > 
		obsIter(radius, m_ObservationImage, m_ObservationImage->GetRequestedRegion());

	NeighborhoodIterator<OutputImageType > 
		confidenceIter(radius, m_ConfidenceImage, m_ConfidenceImage->GetRequestedRegion());

	// setup random iterator
	int totalPixel = 1;
	for( unsigned int dim = 0; dim < OutputImageType::ImageDimension; dim++ )
		totalPixel *= output->GetRequestedRegion().GetSize()[dim];
	ImageRandomConstIteratorWithIndex<OutputImageType> helpIterRandom(output, output->GetRequestedRegion());
	helpIterRandom.SetNumberOfSamples(totalPixel);

	// setup linear iterator
	ImageRegionIterator<OutputImageType> helpIterLinear(output, output->GetRequestedRegion());

	// shadow image
	ImageRegionIterator<OutputImageType> shadow(m_Shadow, m_Shadow->GetRequestedRegion());
	ImageRegionIterator<OutputImageType> outputIter(output, output->GetRequestedRegion());
	if(!m_DuringIterationUpdate) {
		shadow.GoToBegin(); outputIter.GoToBegin();
		for(shadow.GoToBegin(), outputIter.GoToBegin();!outputIter.IsAtEnd(); ++shadow, ++outputIter)
			shadow.Set(outputIter.Get());
	}

	// try new center values
	const EnergyValueType bigPrecalc = 1;
	const EnergyValueType smallPrecalc = 0.1;
	
	VectorImageType::PixelType value0;
	value0[0] = 0; value0[1] = 0;
	VectorImageType::PixelType value1;
	value1[0] = smallPrecalc; value1[1] = 0;
	VectorImageType::PixelType value2;
	value2[0] = -smallPrecalc; value2[1] = 0;
	VectorImageType::PixelType value3;
	value3[0] = bigPrecalc; value3[1] = 0;
	VectorImageType::PixelType value4;
	value4[0] = -bigPrecalc; value4[1] = 0;
	VectorImageType::PixelType value5;
	value5[0] = 0; value5[1] = smallPrecalc;
	VectorImageType::PixelType value6;
	value6[0] = 0; value6[1] = -smallPrecalc;
	VectorImageType::PixelType value7;
	value7[0] = 0; value7[1] = bigPrecalc;
	VectorImageType::PixelType value8;
	value8[0] = 0; value8[1] = -bigPrecalc;

	typedef itk::MinimumDecisionRule DecisionRuleType;
	DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();
	std::vector< double > energyValues;
	VectorImageType::PixelType newCenterPixelChange[9] = { value0, value1, value2,
		value3, value4, value5, value6, value7, value8};

	// loop over neighborhoods and write to output image
	for(helpIterLinear.GoToBegin(); !helpIterLinear.IsAtEnd(); ++helpIterLinear) {
		// set new position of neighborhood iterator
		OutputImageType::IndexType pos = helpIterLinear.GetIndex();
		dispIter.SetLocation(pos);
		labelIter.SetLocation(pos);
		obsIter.SetLocation(pos);
		confidenceIter.SetLocation(pos);

		energyValues.clear();

		for(int j = 0; j < 9; j++) {
			VectorImageType::PixelType newCenterPixel = 
				dispIter.GetCenterPixel() + newCenterPixelChange[j];
			
			// own itk energy function implementation
			float localEnergy;
			float localPriorEnergy;
			float localObservationEnergy;
			(*m_EnergyFunction)(newCenterPixel, dispIter, labelIter, obsIter, confidenceIter, 
				localEnergy, localPriorEnergy, localObservationEnergy);

			energyValues.push_back( localEnergy );

			/*
			float ratio = runningEnergy/localEnergy;
			bool takeIt = false;
			if(ratio >= 1) takeIt = true;
			else {
				float propability = (rand() % 100)/100.0;
				if(ratio > propability) takeIt = true;
			}
			if(takeIt) {*/
			/*
			if(runningEnergy > localEnergy) {
				if(m_DuringIterationUpdate)
					dispIter.SetCenterPixel( newCenterPixel );
				else
					m_Shadow->SetPixel(pos, newCenterPixel);
				runningEnergy = localEnergy;
			}*/
		}
		// take minimum values
		unsigned int lowestEnergy = decisionRule->Evaluate( energyValues );
		VectorImageType::PixelType newValue = dispIter.GetCenterPixel() + newCenterPixelChange[lowestEnergy];

		if(m_DuringIterationUpdate)
			dispIter.SetCenterPixel(newValue);
		else
			m_Shadow->SetPixel(pos, newValue);
	}

	if(!m_DuringIterationUpdate) {
		// now write to output image
		for(shadow.GoToBegin(), outputIter.GoToBegin(); !shadow.IsAtEnd(); ++shadow, ++outputIter)
			outputIter.Set(shadow.Get());
	}

#ifdef _DEBUG
	// just debugging
	if(totalPixel == 1 || totalPixel == 4) {
		for(outputIter.GoToBegin(), labelIter.GoToBegin(); !outputIter.IsAtEnd(); ++outputIter, ++labelIter) {
			std::cout << "displacement " << outputIter.GetIndex() << " " << outputIter.Get() << std::endl;
			energyFile << "displacement " << outputIter.GetIndex() << " " << outputIter.Get() << std::endl;
			std::cout << "label " << labelIter.GetIndex() << " " << labelIter.GetCenterPixel() << std::endl;
			energyFile << "label " << labelIter.GetIndex() << " " << labelIter.GetCenterPixel() << std::endl;
		}
	}
#endif
}

template<class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
float
MRFRegularizationFilter<TInputImage,TOutputImage,TNLabelImage,TEnergyFunction>
::ComputeRMSDFieldEnergy(typename OutputImageType::Pointer output, typename OutputImageType::Pointer backUp,
					 typename NeighborhoodIterator<OutputImageType>::RadiusType radius)
{
	NeighborhoodIterator<OutputImageType > 
		dispIter(radius, output, output->GetRequestedRegion());

	NeighborhoodIterator<OutputImageType >
		backUpIter(radius, backUp, backUp->GetRequestedRegion());

	NeighborhoodIterator<TNLabelImage > 
		labelIter(radius, m_Labels, m_Labels->GetRequestedRegion());

	NeighborhoodIterator<OutputImageType > 
		obsIter(radius, m_ObservationImage, m_ObservationImage->GetRequestedRegion());

	NeighborhoodIterator<OutputImageType > 
		confidenceIter(radius, m_ConfidenceImage, m_ConfidenceImage->GetRequestedRegion());

	// loop over neighborhoods and write to output image
	float rmsd = 0;
	float diff;
	unsigned int n = 0;
	float localEnergyCurrent;
	float localEnergyBackUp;
	float localPriorEnergy;
	float localObservationEnergy;
	for(dispIter.GoToBegin(), backUpIter.GoToBegin(), labelIter.GoToBegin(), obsIter.GoToBegin(), 
			confidenceIter.GoToBegin(); !dispIter.IsAtEnd(); ++dispIter, ++labelIter, ++obsIter, 
			++confidenceIter, ++backUpIter, ++n) {
		// own itk energy function
		(*m_EnergyFunction)(dispIter.GetCenterPixel(), dispIter, labelIter, obsIter, confidenceIter, 
			localEnergyCurrent, localPriorEnergy, localObservationEnergy);
		(*m_EnergyFunction)(backUpIter.GetCenterPixel(), backUpIter, labelIter, obsIter, confidenceIter, 
			localEnergyBackUp, localPriorEnergy, localObservationEnergy);
		diff = localEnergyCurrent - localEnergyBackUp;
		rmsd += diff*diff;
	}

	return std::sqrt(rmsd / n);
}

template<class TInputImage, class TOutputImage, class TNLabelImage, class TEnergyFunction>
float
MRFRegularizationFilter<TInputImage,TOutputImage,TNLabelImage,TEnergyFunction>
::ComputeRMSDDisplacement(typename ImageRegionIterator<OutputImageType> &backUpIter, 
					 typename ImageRegionIterator<OutputImageType> &currentIter)
{
	OutputImageType::PixelType v;
	float rmsd = 0;
	float diff;
	unsigned int n = 0;
	for(backUpIter.GoToBegin(), currentIter.GoToBegin(); !backUpIter.IsAtEnd(); ++backUpIter, 
			++currentIter, ++n) {
		v = backUpIter.Get() - currentIter.Get();
		diff = std::sqrt(v[0]*v[0] + v[1]*v[1]);
		rmsd += diff*diff;
	}

	return std::sqrt(rmsd / n);
}

} // end namespace itk

#endif
