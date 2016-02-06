/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	
 */

#include "Definitions.h"
#include "MetropolisSamplerMinimizer.h"
#include "Utils.h"
#include "MechanicalProperties.h"

#include <iostream>
#include <fstream>
#include <string.h>
#include <windows.h>

#include <vtkAbaqusPointsReader.h>

#include <itkMRFRegularizationFilter.h>
#include <itkYoungPoissonEnergyFunction.h>

#include <itkAmoebaOptimizer.h>
#include <RegularizationCostFunction.h>

using namespace dvfRegularization;

int main(int argc,char *argv[])
{
	// input image vtkImageData's
	vtkImageData* inputImage = vtkImageData::New();
	int dim[3] = {301,301,1};
	//inputImage->SetExtent( 0, width, 0, height, 0, 0);
	inputImage->SetDimensions(dim);
	inputImage->SetOrigin( 0.0, 0.0, 0.0);
	inputImage->SetSpacing( 1.0, 1.0, 1.0);
	inputImage->SetScalarTypeToFloat();
	inputImage->SetNumberOfScalarComponents(1);
	inputImage->AllocateScalars();
	for(int i = 0; i < dim[0]; i++)
	{
		for(int j = 0; j < dim[1]; j++)
		{
			// get scalar pointer to current pixel
			float *currentPixel = (float*) inputImage->GetScalarPointer( i , j , 0 );
			// set scalar value accordingly
			if(j >= 50 && i >= 0 && j <= 150 && i <= 100) *currentPixel = -20;
			else if(j > 50 && i >= 0 && j <= 250 && i <= 200) *currentPixel = 20;
			else *currentPixel = -100;
		}
	}

	// create label
	ImageType::Pointer label = Utils::createLabels(dim);

	// create the corresponding mechanical properties
	vector<MechanicalProperties*> mechanicalPropertiesVector;
	mechanicalPropertiesVector.push_back(new MechanicalProperties(15000000, 0.3));
	mechanicalPropertiesVector.push_back(new MechanicalProperties(10, 0.3));
	mechanicalPropertiesVector.push_back(new MechanicalProperties(20000000, 0.3));

	// create deformation
	ImageType::Pointer elastic[3];
	elastic[0] = Utils::oneDirectionDVF(dim, 30);
	elastic[1] = Utils::twoDirectionDVF(dim, 5);
	elastic[2] = Utils::twoDirectionDVF(dim, 0);
	
	// combine all directions
	VectorImageType::Pointer elasticXYZ = Utils::createVectorImage(elastic);

	// read abaqus results for comparison
	/*
	vtkAbaqusPointsReader* abaqusResults = vtkAbaqusPointsReader::New();
	abaqusResults->SetFileName("../DVFRegularization/label_u1_u2.txt");
	//Utils::visualizeDvfArrow(abaqusResults->GetOutput());
	*/

	// convert image to itk
	/*
	typedef itk::VTKImageToImageFilter<VectorImageType> ITKExporterVectorType;
	ITKExporterVectorType::Pointer ITKExporter = ITKExporterVectorType::New();
	ITKExporter->SetInput(abaqusResults->GetOutput());
	ITKExporter->Update();
	*/

	VectorImageType::Pointer abaqusItkImage = 
	Utils::readAbaqusFile("../DVFRegularization/two_tissue_abaqus.rpy");

	vector<ImageType::Pointer> elasticVector;
	for(unsigned int i = 0; i < 3; i++) elasticVector.push_back(elastic[i]);

	// markov chain (= no memory) monte carlo (= ramdon) method
	MetropolisSamplerMinimizer metropolis;
	metropolis.run(elasticXYZ, elasticVector, label, mechanicalPropertiesVector, abaqusItkImage);

	/*
	// simplex optimization method
	typedef itk::AmoebaOptimizer OptimizerType;
	typedef OptimizerType::InternalOptimizerType  vnlOptimizerType;
	OptimizerType::Pointer itkOptimizer = OptimizerType::New();
	RegularizationCostFunction::Pointer costFunction = RegularizationCostFunction::New();
	costFunction->SetInitial(elasticXYZ);
	costFunction->SetElasticVector(&elasticVector);
	costFunction->SetLabel(label);
	costFunction->SetMechanicalProperties(&mechanicalPropertiesVector);
	costFunction->SetAbaqus(abaqusItkImage);
	itkOptimizer->SetCostFunction(costFunction);
	OptimizerType::ParametersType initialValue(4);
	initialValue[0] = 0.5; initialValue[1] = 1.0; initialValue[2] = 1.0; initialValue[3] = 1.0;
	itkOptimizer->SetInitialPosition( initialValue );
	itkOptimizer->SetParametersConvergenceTolerance( 0.001 );
	itkOptimizer->SetFunctionConvergenceTolerance(0.001); // 0.1%
	itkOptimizer->SetMaximumNumberOfIterations( 200 );

	vnlOptimizerType * vnlOptimizer = itkOptimizer->GetOptimizer();
	try  {
		vnlOptimizer->verbose = true;
		std::cout << "Run for " << itkOptimizer->GetMaximumNumberOfIterations();
		std::cout << " iterations." << std::endl;
		itkOptimizer->StartOptimization();
	}
	catch( itk::ExceptionObject & e ) {
		std::cout << "Exception thrown ! " << std::endl;
		std::cout << "An error ocurred during Optimization" << std::endl;
		std::cout << "Location    = " << e.GetLocation()    << std::endl;
		std::cout << "Description = " << e.GetDescription() << std::endl;
	}
	std::cout << "Number of evals = " << vnlOptimizer->get_num_evaluations() << std::endl;
	std::cout << "Optimizer: " << itkOptimizer;
	*/

	inputImage->Delete();
}
