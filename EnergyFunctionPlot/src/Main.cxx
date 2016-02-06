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
#include "Utils.h"
#include "MechanicalProperties.h"
#include "DisplacementBoundaryCondition.h"

#include <itkYoungPoissonEnergyFunction.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <windows.h>

using namespace dvfRegularization;

int main(int argc,char *argv[])
{
	// creat deformation field image
	VectorImageType::Pointer elasticXYZ = VectorImageType::New();
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0] = 2; size[1] = 1;
	start[0] = 0; start[1] = 0;
	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	elasticXYZ->SetRegions( region );
	elasticXYZ->Allocate();
	// set image to zero
	VectorImageType::IndexType index1, index2;
	VectorImageType::PixelType value1, value2;
	index1[0] = 0; index1[1] = 0;
	index2[0] = 1; index2[1] = 0;
	value1[0] = 0; value1[1] = 0;
	value2[0] = 0; value2[1] = 0;
	elasticXYZ->SetPixel(index1, value1);
	elasticXYZ->SetPixel(index2, value2);

	// create label image
	FloatImageType::Pointer labels = FloatImageType::New();
	labels->SetRegions(region);
	labels->Allocate();
	labels->SetPixel(index1, 28.9);
	labels->SetPixel(index2, 13.4);

	// neighborhood radius
	itk::NeighborhoodIterator<VectorImageType>::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;

	typedef DisplacementBoundaryCondition<VectorImageType> BC;
	BC bc;
	VectorPixelType constant;
	constant[0] = 0; constant[1] = 0;
	VectorPixelType displacementRight;
	displacementRight[0] = -15; displacementRight[1] = 0;
	VectorPixelType displacementLeft;
	displacementLeft[0] = 0; displacementLeft[1] = 0;
	VectorPixelType displacementTop;
	displacementTop[0] = 0; displacementTop[1] = 0;
	VectorPixelType displacementBottom;
	displacementBottom[0] = 0; displacementBottom[1] = 0;
	bc.SetConstant(constant);
	bc.SetMaxDisplacementRight(displacementRight);
	bc.SetMaxDisplacementLeft(displacementLeft);
	bc.SetMaxDisplacementTop(displacementTop);
	bc.SetMaxDisplacementBottom(displacementBottom);

	NeighborhoodIterator<VectorImageType,BC> 
		dispIter(radius, elasticXYZ, elasticXYZ->GetRequestedRegion());
	dispIter.OverrideBoundaryCondition(&bc);

	NeighborhoodIterator<FloatImageType,ZeroFluxNeumannBoundaryCondition<FloatImageType> > 
		labelIter(radius, labels, labels->GetRequestedRegion());

	itk::YoungPoissonEnergyFunction<VectorImageType,FloatImageType,BC> energyFunction(1.0, 1.0, 1.0);

	ofstream xyPlotFile;
	xyPlotFile.open("x_y_energy_value.txt");

	value1[0] = 0; value1[1] = 0;
	value2[0] = 0; value2[1] = 0;
	for(int x1 = 0; x1 <= 30; x1 += 1) {
		value1[0] = -x1;
		elasticXYZ->SetPixel(index1, value1);
		for(int x2 = 0; x2 <= 30; x2 += 1) {
			value2[0] = -x2;
			elasticXYZ->SetPixel(index2, value2);
			float totalEnergy = 0;
			for(dispIter.GoToBegin(), labelIter.GoToBegin(); !dispIter.IsAtEnd(); ++dispIter, ++labelIter) {
				float localEnergy;
				float localPriorEnergy;
				float localObservationEnergy;
				energyFunction(dispIter.GetCenterPixel(), dispIter, labelIter, 
					localEnergy, localPriorEnergy, localObservationEnergy);
				totalEnergy += localEnergy;
			}
			xyPlotFile << value1[0] << " " << value2[0] << " " << totalEnergy << std::endl;
		}
	}

	xyPlotFile.close();
}
