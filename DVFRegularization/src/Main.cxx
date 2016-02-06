/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Represents an energy function with its two parts.
 *			Umodel + Uprior = Utotal
 */

#include "Definitions.h"
#include "Utils.h"
#include "MechanicalProperties.h"
#include "DisplacementBoundaryCondition.h"

#include "itkMRFRegularizationFilter.h"
#include "itkMRFRegularizationFilterGPU.h"
#include "itkYoungPoissonEnergyFunction.h"
#include "itkInterpolateNextStepFilter.h"

#include <itkImageRegionIterator.h>
#include <itkWarpImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageFileReader.h>
#include <itkResampleImageFilter.h>
#include <itkFlipImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkPointSet.h>
#include <itkPointSetToImageFilter.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <windows.h>

using namespace dvfRegularization;

int main(int argc,char *argv[])
{
	if(argc < 5) {
		std::cout<<"USAGE: tissueImage.png labelImage.png bc.mhd ATI=0|NVIDIA=1 abaqus.mhd"<<std::endl;
		Sleep(5000);
		exit(-1);
	}

	// image file reader
	typedef itk::ImageFileReader< Vector3ImageType > Reader3DType;
	typedef itk::ImageFileWriter< Vector3ImageType > Writer3DType;
	Writer3DType::Pointer writer3d = Writer3DType::New();
	typedef itk::ImageFileWriter< GreyImageType > Writer2DType;
	Writer2DType::Pointer writer2d = Writer2DType::New();
	typedef itk::ImageFileReader<GreyImageType> ImageReaderType;
	typedef itk::ImageFileWriter< FloatImageType > Writer1DType;
	Writer1DType::Pointer writer1d = Writer1DType::New();
	typedef itk::ImageFileReader<GreyImageType> ImageReaderType;

	// biomechanical properties
	ImageReaderType::Pointer tissueReader = ImageReaderType::New();
	tissueReader->SetFileName(argv[1]);
	tissueReader->Update();
	GreyImageType::Pointer tissueImage = tissueReader->GetOutput();
	ImageReaderType::Pointer labelReader = ImageReaderType::New();
	labelReader->SetFileName(argv[2]);
	labelReader->Update();
	GreyImageType::Pointer labelImage = labelReader->GetOutput();
	VectorImageType::Pointer elasticXYZ;
	vector<MechanicalProperties*> mechanicalPropertiesVector;

	// 3d dvf for paraview
	Vector3ImageType::Pointer paraview = Vector3ImageType::New();
	VectorImageType::RegionType  region2D = tissueImage->GetBufferedRegion();
	VectorImageType::IndexType   index2D  = region2D.GetIndex();
	VectorImageType::SizeType    size2D   = region2D.GetSize();
	Vector3ImageType::RegionType  region3D;
	Vector3ImageType::IndexType   index3D;
	Vector3ImageType::SizeType    size3D;
	index3D[0] = index2D[0];
	index3D[1] = index2D[1];
	index3D[2] = 0;
	size3D[0]  = size2D[0];
	size3D[1]  = size2D[1];
	size3D[2]  = 1;
	region3D.SetSize( size3D );
	region3D.SetIndex( index3D );
	paraview->SetRegions( region3D );
	paraview->Allocate();
	itk::ImageRegionIterator<Vector3ImageType> iterParaview(paraview, paraview->GetRequestedRegion());

	// abaqus file name for comparison
	std::string abaqusFile;

	// create the corresponding mechanical properties, order is important
	// second muslce tissue
	mechanicalPropertiesVector.push_back(new MechanicalProperties(10, 0.3));
	// muscle tissue
	mechanicalPropertiesVector.push_back(new MechanicalProperties(20, 0.3));
	// bony tissue
	mechanicalPropertiesVector.push_back(new MechanicalProperties(40, 0.3));

	// do hierachy optimization
	int pixelCount = 1;
	for(int i = 0; i < tissueImage->GetRequestedRegion().GetImageDimension(); ++i)
		pixelCount *= tissueImage->GetRequestedRegion().GetSize()[i];
	int levels = std::log((float)pixelCount)/std::log((float)4) + 2;
	// definitions
	typedef itk::MRFRegularizationFilterGPU<VectorImageType,VectorImageType,FloatImageType> MRFFilterGPU;
	MRFFilterGPU::Pointer mrfFilterGPU;
	typedef itk::MRFRegularizationFilter<VectorImageType,VectorImageType,FloatImageType,
		itk::YoungPoissonEnergyFunction<VectorImageType,FloatImageType> > MRFFilter;
	MRFFilter::Pointer mrfFilter;
	
	typedef itk::NearestNeighborInterpolateImageFunction<FloatImageType, double > InterpolatorType;
	//typedef itk::LinearInterpolateImageFunction<FloatImageType, double > InterpolatorType;
	typedef InterpolateNextStepFilter<VectorImageType,VectorImageType,InterpolatorType> NextStepFilter;
	NextStepFilter::Pointer nextStepFilter;
	NextStepFilter::Pointer confidenceInterpolation;

	typedef itk::LinearInterpolateImageFunction<FloatImageType, double > ElasticInterpolatorType;
	typedef InterpolateNextStepFilter<VectorImageType,VectorImageType,ElasticInterpolatorType> 
		ElasticInterpolationFilter;
	ElasticInterpolationFilter::Pointer elasticInterpolationFilter;
	
	// boundary conditions as observation image
	Reader3DType::Pointer reader3d = Reader3DType::New();
	reader3d->SetFileName(argv[3]);
	reader3d->Update();
	VectorImageType::Pointer observationImage = VectorImageType::New();
	observationImage->SetRegions(tissueImage->GetRequestedRegion());
	observationImage->Allocate();
	itk::ImageRegionIterator<VectorImageType> obsIter(observationImage, 
		observationImage->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> bcIter(reader3d->GetOutput(), 
		reader3d->GetOutput()->GetRequestedRegion());
	for(obsIter.GoToBegin(), bcIter.GoToBegin(); !obsIter.IsAtEnd(); ++obsIter, ++bcIter) {
		VectorImageType::PixelType value;
		value[0] = bcIter.Get()[0]; value[1] = bcIter.Get()[1];
		obsIter.Set(value);
	}
	// confidence image
	Reader3DType::Pointer confidenceReader = Reader3DType::New();
	confidenceReader->SetFileName(argv[4]);
	confidenceReader->Update();
	VectorImageType::Pointer confidenceImage = VectorImageType::New();
	confidenceImage->SetRegions(tissueImage->GetRequestedRegion());
	confidenceImage->Allocate();
	itk::ImageRegionIterator<VectorImageType> confidenceIter(confidenceImage, 
		confidenceImage->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> readerIter(confidenceReader->GetOutput(), 
		confidenceReader->GetOutput()->GetRequestedRegion());
	for(confidenceIter.GoToBegin(), readerIter.GoToBegin(); !confidenceIter.IsAtEnd(); 
			++confidenceIter, ++readerIter) {
		VectorImageType::PixelType value;
		value[0] = readerIter.Get()[0]; value[1] = readerIter.Get()[1];
		confidenceIter.Set(value);
	}

	// energy function
	itk::YoungPoissonEnergyFunction<VectorImageType,FloatImageType> energyFunction(0.5, 1.0, 1.0);

	const bool hierachyEnabled = true;
	const bool gpu = true;

	if(gpu) {
		// workaround
		initGPU();
	}

	// timing vars 
	double start, end;
	// calc on CPU
	start=clock();

	for(int level = 0; level < levels; ++level) {
		// prepare level
#ifdef _DEBUG
		std::cout << "--- Start with level " << level << " ---" << std::endl;
#endif
		int currentPixelCount;
		if(pixelCount > std::pow((float)4,level))
			currentPixelCount = std::pow((float)4,level);
		else
			currentPixelCount = pixelCount;

		if(gpu) {
			mrfFilterGPU = MRFFilterGPU::New();
			// "ati=0" or "nvidia=1"
			mrfFilterGPU->SetGPUType(atoi(argv[5]));
			mrfFilterGPU->SetNumberOfIterations(100);
		}
		else {
			// prepare markov filter
			mrfFilter = MRFFilter::New();
			mrfFilter->SetMaximumNumberOfIterations(100);
			mrfFilter->SetDuringIterationUpdate(false);
			mrfFilter->SetStopCriteria(0.0); //0.0001
			mrfFilter->SetEnergyFunction(&energyFunction);
			mrfFilter->SetVisualize(false);
			mrfFilter->SetRandomIteration(true);
		}

		// interpolate initial displacement vector field based on last iteration
		if(hierachyEnabled) {
			if(gpu) {
				mrfFilterGPU->SetLabels(Utils::extractYM(Utils::combineLabels(currentPixelCount, labelImage),
					mechanicalPropertiesVector));
			}
			else {
				mrfFilter->SetLabels(Utils::extractYM(Utils::combineLabels(currentPixelCount, labelImage),
					mechanicalPropertiesVector));
			}
			//elasticXYZ = Utils::extrapolateNearestNeighborInitialDVF(currentPixelCount, elasticXYZ);
			if(currentPixelCount == 1) {
				elasticXYZ = VectorImageType::New();
				VectorImageType::IndexType start;
				VectorImageType::SizeType  size;
				size[0] = 1; size[1] = 1;
				start[0] = 0; start[1] = 0;
				VectorImageType::RegionType region;
				region.SetSize( size );
				region.SetIndex( start );
				elasticXYZ->SetRegions( region );
				elasticXYZ->Allocate();
				VectorImageType::PixelType value;
				value[0] = 0; value[1] = 0;
				elasticXYZ->FillBuffer(value);
			}
			else {
				elasticInterpolationFilter = ElasticInterpolationFilter::New();
				elasticInterpolationFilter->SetNumberOfPixels(currentPixelCount);
				elasticInterpolationFilter->SetInput(elasticXYZ);
				elasticInterpolationFilter->Update();
				elasticXYZ = elasticInterpolationFilter->GetOutput();
			}
			nextStepFilter = NextStepFilter::New();
			nextStepFilter->SetNumberOfPixels(currentPixelCount);
			nextStepFilter->SetInput(observationImage);
			nextStepFilter->Update();
			confidenceInterpolation = NextStepFilter::New();
			confidenceInterpolation->SetNumberOfPixels(currentPixelCount);
			confidenceInterpolation->SetInput(confidenceImage);
			confidenceInterpolation->Update();
		}
		else {
			if(gpu) {
				mrfFilterGPU->SetLabels(Utils::extractYM(Utils::combineLabels(pixelCount, labelImage),
					mechanicalPropertiesVector));
			}
			else {
				mrfFilter->SetLabels(Utils::extractYM(Utils::combineLabels(pixelCount, labelImage),
					mechanicalPropertiesVector));
			}
			elasticXYZ = VectorImageType::New();
			VectorImageType::IndexType start;
			VectorImageType::SizeType  size;
			size[0] = std::sqrt((float)pixelCount); size[1] = std::sqrt((float)pixelCount);
			start[0] = 0; start[1] = 0;
			VectorImageType::RegionType region;
			region.SetSize( size );
			region.SetIndex( start );
			elasticXYZ->SetRegions( region );
			elasticXYZ->Allocate();
			VectorImageType::PixelType value;
			value[0] = 0; value[1] = 0;
			elasticXYZ->FillBuffer(value);
			// leave for loop
			level = levels;
		}

		if(gpu) {
			if(hierachyEnabled) {
				mrfFilterGPU->SetObservationImage(nextStepFilter->GetOutput());
				mrfFilterGPU->SetConfidenceImage(confidenceInterpolation->GetOutput());
			}
			else {
				mrfFilterGPU->SetObservationImage(observationImage);
				mrfFilterGPU->SetConfidenceImage(confidenceImage);
			}
			mrfFilterGPU->SetInput(elasticXYZ);
			mrfFilterGPU->Update();
		}
		else {
			if(hierachyEnabled) {
                mrfFilter->SetObservationImage(nextStepFilter->GetOutput());
				mrfFilter->SetConfidenceImage(confidenceInterpolation->GetOutput());
			}
			else {
				mrfFilter->SetObservationImage(observationImage);
				mrfFilter->SetConfidenceImage(confidenceImage);
			}
			mrfFilter->SetInput(elasticXYZ);
			mrfFilter->Update();
		}

		// copy
		if(gpu) {
			itk::ImageRegionIterator<VectorImageType> iterFrom(mrfFilterGPU->GetOutput(), 
				mrfFilterGPU->GetOutput()->GetRequestedRegion());
			itk::ImageRegionIterator<VectorImageType> iterTo(elasticXYZ, elasticXYZ->GetRequestedRegion());
			for(iterFrom.GoToBegin(), iterTo.GoToBegin(); !iterFrom.IsAtEnd(); ++iterFrom, ++iterTo)
				iterTo.Set(iterFrom.Get());
		}
		else {
			itk::ImageRegionIterator<VectorImageType> iterFrom(mrfFilter->GetOutput(), 
				mrfFilter->GetOutput()->GetRequestedRegion());
			itk::ImageRegionIterator<VectorImageType> iterTo(elasticXYZ, elasticXYZ->GetRequestedRegion());
			for(iterFrom.GoToBegin(), iterTo.GoToBegin(); !iterFrom.IsAtEnd(); ++iterFrom, ++iterTo)
				iterTo.Set(iterFrom.Get());
		}

#ifdef _DEBUG

		// matlab analysis
		itk::ImageRegionIterator<VectorImageType> dispIter(elasticXYZ, elasticXYZ->GetRequestedRegion());
		
		// save data horizontal line to text file for visualizing with matlab
		std::stringstream hPlotName;
		hPlotName << "h_plot_" << level << ".txt";
		ofstream hPlotFile;
		hPlotFile.open(hPlotName.str().c_str());
		for(dispIter.GoToBegin(); !dispIter.IsAtEnd(); ++dispIter) {
			//if(dispIter.GetIndex()[1] == 5*elasticXYZ->GetRequestedRegion().GetSize()[1]/6) {
			if(dispIter.GetIndex()[1] == elasticXYZ->GetRequestedRegion().GetSize()[1]/2) {
				hPlotFile << dispIter.GetIndex()[0] << " " << dispIter.Get()[0] <<
					" " << dispIter.Get()[1] << std::endl;
			}
		}
		hPlotFile.close();

		// save data vertical line to text file for visualizing with matlab
		std::stringstream vPlotName;
		vPlotName << "v_plot_" << level << ".txt";
		ofstream vPlotFile;
		vPlotFile.open(vPlotName.str().c_str());
		// flip y index because itk vtk difference
		unsigned int ySize = elasticXYZ->GetRequestedRegion().GetSize()[1]-1;
		for(dispIter.GoToBegin(); !dispIter.IsAtEnd(); ++dispIter) {
			//if(dispIter.GetIndex()[0] == elasticXYZ->GetRequestedRegion().GetSize()[0]/6) {
			if(dispIter.GetIndex()[0] == elasticXYZ->GetRequestedRegion().GetSize()[0]/2) {
				vPlotFile << ySize-dispIter.GetIndex()[1] << " " << dispIter.Get()[0] <<
					" " << dispIter.Get()[1] << std::endl;
			}
		}
		vPlotFile.close();

		// save data vertical line to text file for visualizing with matlab
		std::stringstream d1PlotName;
		d1PlotName << "d1_plot_" << level << ".txt";
		ofstream d1PlotFile;
		d1PlotFile.open(d1PlotName.str().c_str());
		// flip y index because itk vtk difference
		for(dispIter.GoToBegin(); !dispIter.IsAtEnd(); ++dispIter) {
			if(dispIter.GetIndex()[0] == dispIter.GetIndex()[1]) {
				d1PlotFile << dispIter.GetIndex()[0] << " " << dispIter.Get()[0] <<
					" " << dispIter.Get()[1] << std::endl;
			}
		}
		d1PlotFile.close();

		// save data vertical line to text file for visualizing with matlab
		std::stringstream d2PlotName;
		d2PlotName << "d2_plot_" << level << ".txt";
		ofstream d2PlotFile;
		d2PlotFile.open(d2PlotName.str().c_str());
		// flip y index because itk vtk difference
		for(dispIter.GoToBegin(); !dispIter.IsAtEnd(); ++dispIter) {
			if(dispIter.GetIndex()[0] == ySize-dispIter.GetIndex()[1]) {
				d2PlotFile << dispIter.GetIndex()[0] << " " << dispIter.Get()[0] <<
					" " << dispIter.Get()[1] << std::endl;
			}
		}
		d2PlotFile.close();

#endif

	}

	end = clock();
	double total = (end-start)/CLOCKS_PER_SEC;
	std::cout << "Computation time: " << total << " sec." << std::endl;

	if(gpu) {
		// workaround
		cleanUpGPU();
	}

	// compute total field energy
	NeighborhoodIterator<VectorImageType>::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	NeighborhoodIterator<VectorImageType> dispIter(radius, elasticXYZ, elasticXYZ->GetRequestedRegion());
	FloatImageType::Pointer tmpImage = Utils::extractYM(
		Utils::combineLabels(pixelCount, labelImage), mechanicalPropertiesVector);
	NeighborhoodIterator<FloatImageType> labelIter(radius, tmpImage, tmpImage->GetRequestedRegion());
	NeighborhoodIterator<VectorImageType> neighbObsIter(radius, observationImage, 
		observationImage->GetRequestedRegion());
	NeighborhoodIterator<VectorImageType> neighConfidenceIter(radius, confidenceImage, 
		confidenceImage->GetRequestedRegion());
	// loop over neighborhoods and write to output image
	float meanFieldEnergy = 0;
	float localEnergy = 0;
	unsigned int n = 0;
	float localPriorEnergy;
	float localObservationEnergy;
	for(dispIter.GoToBegin(), labelIter.GoToBegin(), neighbObsIter.GoToBegin(), 
			neighConfidenceIter.GoToBegin(); !dispIter.IsAtEnd(); ++dispIter, ++labelIter, 
			++neighbObsIter, ++neighConfidenceIter, ++n) {
		// own itk energy function
		energyFunction(dispIter.GetCenterPixel(), dispIter, labelIter, neighbObsIter, neighConfidenceIter,
			localEnergy, localPriorEnergy, localObservationEnergy);
		meanFieldEnergy += localEnergy;
	}
	meanFieldEnergy = meanFieldEnergy / n;
	std::cout << "Mean field energy: " << meanFieldEnergy << std::endl;

	// write displacement field to file for further analysis in paraview
	// vtk file extension: .vti
	// meta file extension: .mhd or .mha
	// difference between fem and markov
	itk::ImageRegionIterator<VectorImageType> iterDVF(elasticXYZ, 
		elasticXYZ->GetRequestedRegion());
	for(iterDVF.GoToBegin(), iterParaview.GoToBegin(); !iterDVF.IsAtEnd(); ++iterDVF, ++iterParaview) {
		Vector3ImageType::PixelType vector3d;
		VectorImageType::PixelType vector2d = iterDVF.Get();
		vector3d[0] = vector2d[0];
		vector3d[1] = vector2d[1];
		vector3d[2] = 0;
		iterParaview.Set(vector3d);
	}
	// flip for vtk
	typedef itk::FlipImageFilter<Vector3ImageType> FlipFilter;
	FlipFilter::Pointer flipFilter = FlipFilter::New();
	FlipFilter::FlipAxesArrayType axis;
	axis[0] = 0; axis[1] = 1; axis[2] = 0;
	flipFilter->SetFlipAxes(axis);
	flipFilter->SetInput(paraview);
	writer3d->SetInput(flipFilter->GetOutput());
	writer3d->SetFileName( "displacement_field.mhd" );
    writer3d->Update();

	// compare sampled mrf and abaqus images
	Reader3DType::Pointer abaqusReader = Reader3DType::New();
	abaqusReader->SetFileName(argv[6]);
	abaqusReader->Update();
	// read binary image
	typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;
	BinaryImageReaderType::Pointer binaryReader = BinaryImageReaderType::New();
	binaryReader->SetFileName(argv[7]);
	binaryReader->Update();
	float rmsd;
	Utils::rootMeanSquareDeviation(abaqusReader->GetOutput(), flipFilter->GetOutput(),
		binaryReader->GetOutput(), rmsd);
	//Utils::printDifference(elasticXYZ, mrfFilter->GetOutput());
	cout << "RMSD ABAQUS: " << rmsd << endl;
	typedef itk::SubtractImageFilter<Vector3ImageType,Vector3ImageType,Vector3ImageType> 
		DifferenceFilterType;
	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();
	difference->SetInput1( abaqusReader->GetOutput() );
	difference->SetInput2( flipFilter->GetOutput() );
	difference->Update();
	Vector3ImageType::PixelType value;
	value.Fill(0.0);
	itk::ImageRegionIterator<Vector3ImageType> diffIter(difference->GetOutput(), 
		difference->GetOutput()->GetRequestedRegion());
	itk::ImageRegionIterator<BinaryImageType> binIter(binaryReader->GetOutput(), 
		binaryReader->GetOutput()->GetRequestedRegion());
	for(diffIter.GoToBegin(), binIter.GoToBegin(); !diffIter.IsAtEnd(); ++diffIter, ++binIter) {
		if(binIter.Get() == 0.0)
			diffIter.Set(value);
	}
	writer3d->SetInput(difference->GetOutput());
	writer3d->SetFileName( "fem_mrf_difference.mhd" );
	writer3d->Update();


#ifdef _DEBUG

	// matlab analysis
	itk::ImageRegionIterator<Vector3ImageType> abaqusIter(abaqusReader->GetOutput(), 
		abaqusReader->GetOutput()->GetRequestedRegion());
	
	// save data horizontal line to text file for visualizing with matlab
	ofstream hPlotAFile;
	hPlotAFile.open("h_plot_abaqus.txt");
	for(abaqusIter.GoToBegin(); !abaqusIter.IsAtEnd(); ++abaqusIter) {
		//if(abaqusIter.GetIndex()[1] == 5*abaqusReader->GetOutput()->GetRequestedRegion().GetSize()[1]/6) {
		if(abaqusIter.GetIndex()[1] == abaqusReader->GetOutput()->GetRequestedRegion().GetSize()[1]/2) {
			hPlotAFile << abaqusIter.GetIndex()[0] << " " << abaqusIter.Get()[0] <<
				" " << abaqusIter.Get()[1] << std::endl;
		}
	}
	hPlotAFile.close();

	// save data vertical line to text file for visualizing with matlab
	ofstream vPlotAFile;
	vPlotAFile.open("v_plot_abaqus.txt");
	for(abaqusIter.GoToBegin(); !abaqusIter.IsAtEnd(); ++abaqusIter) {
		//if(abaqusIter.GetIndex()[0] == elasticXYZ->GetRequestedRegion().GetSize()[0]/6) {
		if(abaqusIter.GetIndex()[0] == abaqusReader->GetOutput()->GetRequestedRegion().GetSize()[0]/2) {
			vPlotAFile << abaqusIter.GetIndex()[1] << " " << abaqusIter.Get()[0] <<
				" " << abaqusIter.Get()[1] << std::endl;
		}
	}
	vPlotAFile.close();

	// save data vertical line to text file for visualizing with matlab
	ofstream d1PlotAFile;
	d1PlotAFile.open("d1_plot_abaqus.txt");
	unsigned int ySizeA = abaqusReader->GetOutput()->GetRequestedRegion().GetSize()[1]-1;
	for(abaqusIter.GoToBegin(); !abaqusIter.IsAtEnd(); ++abaqusIter) {
		if(abaqusIter.GetIndex()[0] == ySizeA-abaqusIter.GetIndex()[1]) {
			d1PlotAFile << abaqusIter.GetIndex()[0] << " " << abaqusIter.Get()[0] <<
				" " << abaqusIter.Get()[1] << std::endl;
		}
	}
	d1PlotAFile.close();

	// save data vertical line to text file for visualizing with matlab
	ofstream d2PlotAFile;
	d2PlotAFile.open("d2_plot_abaqus.txt");
	for(abaqusIter.GoToBegin(); !abaqusIter.IsAtEnd(); ++abaqusIter) {
		if(abaqusIter.GetIndex()[0] == abaqusIter.GetIndex()[1]) {
			d2PlotAFile << abaqusIter.GetIndex()[0] << " " << abaqusIter.Get()[0] <<
				" " << abaqusIter.Get()[1] << std::endl;
		}
	}
	d2PlotAFile.close();

#endif

	// apply deformation with itk
	typedef itk::WarpImageFilter<GreyImageType, GreyImageType, VectorImageType > WarpFilterType;
	WarpFilterType::Pointer warpFilter = WarpFilterType::New();
	typedef itk::LinearInterpolateImageFunction<GreyImageType,double> InterpolatorGreyType;
	InterpolatorGreyType::Pointer interpolatorGrey = InterpolatorGreyType::New();
	warpFilter->SetInterpolator(interpolatorGrey);
	warpFilter->SetOutputSpacing(tissueImage->GetSpacing());
	warpFilter->SetOutputOrigin(tissueImage->GetOrigin());
	warpFilter->SetDeformationField(elasticXYZ);
	warpFilter->SetInput(tissueImage);
	writer2d->SetInput(warpFilter->GetOutput());
	writer2d->SetFileName("deformed_tissue_image.png");
    writer2d->Update();
	// deform grid
	ImageReaderType::Pointer gridReader = ImageReaderType::New();
	gridReader->SetFileName(argv[8]);
	warpFilter->SetInput(gridReader->GetOutput());
	writer2d->SetInput(warpFilter->GetOutput());
	writer2d->SetFileName("deformed_grid_image.png");
	writer2d->Update();

	// apply deformation with itk
	flipFilter->SetInput(abaqusReader->GetOutput());
	flipFilter->Update();
	VectorImageType::Pointer abaqusDeformationImage = VectorImageType::New();
	abaqusDeformationImage->SetRegions(tissueImage->GetRequestedRegion());
	abaqusDeformationImage->Allocate();
	itk::ImageRegionIterator<VectorImageType> againIter(abaqusDeformationImage, 
		abaqusDeformationImage->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> again2Iter(flipFilter->GetOutput(), 
		flipFilter->GetOutput()->GetRequestedRegion());
	for(againIter.GoToBegin(), again2Iter.GoToBegin(); !againIter.IsAtEnd(); ++againIter, ++again2Iter) {
		VectorImageType::PixelType value;
		value[0] = again2Iter.Get()[0]; value[1] = again2Iter.Get()[1];
		againIter.Set(value);
	}
	warpFilter->SetDeformationField(abaqusDeformationImage);
	warpFilter->SetInput(tissueImage);
	writer2d->SetInput(warpFilter->GetOutput());
	writer2d->SetFileName("deformed_tissue_image_abaqus.png");
    writer2d->Update();
	// deform grid
	warpFilter->SetInput(gridReader->GetOutput());
	writer2d->SetInput(warpFilter->GetOutput());
	writer2d->SetFileName("deformed_grid_image_abaqus.png");
	writer2d->Update();

	// apply deformation with itk
	Reader3DType::Pointer imageJReader = Reader3DType::New();
	imageJReader->SetFileName("../ImageJToImage/tumor_displacement_imagej.mhd");
	imageJReader->Update();
	VectorImageType::Pointer imageJDeformationImage = VectorImageType::New();
	imageJDeformationImage->SetRegions(tissueImage->GetRequestedRegion());
	imageJDeformationImage->Allocate();
	itk::ImageRegionIterator<VectorImageType> sourceIter(imageJDeformationImage, 
		imageJDeformationImage->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> dvfIter(imageJReader->GetOutput(), 
		imageJReader->GetOutput()->GetRequestedRegion());
	for(sourceIter.GoToBegin(), dvfIter.GoToBegin(); !sourceIter.IsAtEnd(); ++sourceIter, ++dvfIter) {
		VectorImageType::PixelType value;
		value[0] = dvfIter.Get()[0]; value[1] = dvfIter.Get()[1];
		sourceIter.Set(value);
	}
	warpFilter->SetDeformationField(imageJDeformationImage);
	warpFilter->SetInput(tissueImage);
	writer2d->SetInput(warpFilter->GetOutput());
	writer2d->SetFileName("deformed_tissue_image_imagej.png");
    writer2d->Update();
	// deform grid
	warpFilter->SetInput(gridReader->GetOutput());
	writer2d->SetInput(warpFilter->GetOutput());
	writer2d->SetFileName("deformed_grid_image_imagej.png");
	writer2d->Update();
}
