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

#include "MetropolisSamplerMinimizer.h"
#include <itkEnergyFunction.h>
#include <vtkImageData.h>
#include <limits>
#include "itkMRFRegularizationFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkYoungPoissonEnergyFunction.h"
#include "itkConstantBoundaryCondition.h"

using namespace dvfRegularization;

void MetropolisSamplerMinimizer::run(VectorImageType* inital, vector<ImageType::Pointer>& elasticVector,
									 ImageType* label, vector<MechanicalProperties*>& mechanicalProperties,
									 VectorImageType* abaqus) {
	cout << "Start MetropolisSamplerMinimizer" << std::endl;

	const unsigned int iterations = 500;
	// inital parameters
	double T = 1.0;
	float newAlpha = 0.5;
	float oldAlpha = 0.5;
	float beta0 = 1.0;
	float beta1 = 1.0;
	float gamma = 1.0;
	// modify parameters
	float alphaCorrections[4] = {0.001, -0.001, 0.1, -0.1};
	float beta0Corrections[4] = {1.0, -1.0, 10, -10};
	float beta1Corrections[4] = {1.0, -1.0, 10, -10};
	float gammaCorrections[4] = {1.0, -1.0, 10, -10};
	double oldEnergy = numeric_limits<double>::max();
	cout << "Start with maximum energy = " << oldEnergy << endl;

	// write energy value to a file for matlab analysis
	ofstream energyFile;
	energyFile.open ("parameterFile.txt");

	for(unsigned int iteration = 0; iteration < iterations; iteration++) {
		cout << "start iteration = " << iteration << endl;

		cout << "run regularization with the following parameter:" << endl;
		cout << "alpha = " << newAlpha << endl;
		cout << "beta0 = " << beta0 << endl;
		cout << "beta1 = " << beta1 << endl;
		cout << "gamma = " << gamma << endl;
		// run markov own itk filter
		typedef itk::ConstantBoundaryCondition<VectorImageType> BC;
		//typedef itk::ZeroFluxNeumannBoundaryCondition<ImageType> BC;
		typedef itk::MRFRegularizationFilter<VectorImageType,VectorImageType,
			itk::YoungPoissonEnergyFunction<VectorImageType,BC>, BC> MRFFilter;
		MRFFilter::Pointer mrfFilter = MRFFilter::New();
		mrfFilter->SetInput(inital);
		mrfFilter->SetMaximumNumberOfIterations(50);
		mrfFilter->SetDuringIterationUpdate(false);
		mrfFilter->SetStopCriteria(0.01);
		mrfFilter->SetLabels(label);
		mrfFilter->SetMechanicalProperties(&mechanicalProperties);
		ImageType::Pointer tmpA[3];
		for(unsigned int i = 0; i < 3; i++) tmpA[i] = elasticVector[i];
		itk::YoungPoissonEnergyFunction<VectorImageType,BC> energyFunction(newAlpha, beta0, beta1, gamma, 
			tmpA);
		mrfFilter->SetEnergyFunction(&energyFunction);
		mrfFilter->SetVisualize(false);
		mrfFilter->Update();

		double newEnergy = Utils::leastSquare(abaqus, Utils::sampleDVF(mrfFilter->GetOutput()));
		// compare sampled mrf and abaqus images
		cout << "new least square difference = " << newEnergy << endl;

		double delta = newEnergy - oldEnergy;
		cout << "delta = " << delta << endl;
		cout << "exp = " << exp(-delta/T) << endl;
		double P = std::min(1.0,exp(-delta/T));
		if((rand() % 101)/100.0f <= P) {
			cout << "take new energy" << endl;
			oldEnergy = newEnergy;
			// write results to file
			energyFile << iteration << " " << newEnergy << " " << newAlpha << " "
				<< beta0 << " " << beta1 << " " << gamma << endl;
		}
		else {
			cout << "keep old alpha" << endl;
			newAlpha = oldAlpha;
		}
		
		// choose parameter to modify for next iteration
		//unsigned int chooseParameter = rand() % 4;
		unsigned int chooseParameter = 0;
		unsigned int alphaIndex = rand() % 4;
		unsigned int otherIndex = rand() % 4;
		switch(chooseParameter) {
			case 0:
				// alpha
				newAlpha += alphaCorrections[alphaIndex];
				cout << "modify alpha by " << alphaCorrections[alphaIndex] << endl;
				if(newAlpha < 0) newAlpha = 0;
				else if(newAlpha > 1) newAlpha = 1;
				break;
			case 1:
				// beta0
				beta0 += beta0Corrections[otherIndex];
				cout << "modify beta0 by " << beta0Corrections[otherIndex] << endl;
				if(beta0 < 0) beta0 = 0;
				break;
			case 2:
				// beta1
				beta1 += beta1Corrections[otherIndex];
				cout << "modify beta1 by " << beta1Corrections[otherIndex] << endl;
				if(beta1 < 0) beta1 = 0;
				break;
			case 3:
				// gamma
				gamma += gammaCorrections[otherIndex];
				cout << "modify gamma by " << gammaCorrections[otherIndex] << endl;
				if(gamma < 0) gamma = 0;
				break;
		}

		// annealing
		//T = T /log((100.0+iteration))*log(100.0);
		cout << "New temperature = " << T << endl;

		cout << "------------------------------------------------" << endl;
	}
}
