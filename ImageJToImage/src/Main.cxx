/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	
 */

#include "Utils.h"
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkFlipImageFilter.h>
#include "itkInterpolateNextStepFilter.h"
#include <itkPointSet.h>
#include <itkPointSetToImageFilter.h>
#include <itkMatrix.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

int main(int argc,char *argv[])
{
	// dvf
	typedef itk::Vector<float,3> Image3dPixelType;
	typedef itk::Image<Image3dPixelType,3> Image3dType;

	// image dimension
	const int xDim = 300;
	const int yDim = 300;
	const int zDim = 1;

	// save to file
	typedef itk::ImageFileWriter< Image3dType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	// load fem results with rectangular mesh

	std::ifstream rawFile("raw_imagej_dvf.txt");
	if(!rawFile) {
		cerr << "raw_imagej_dvf.txt file not found" << endl;
		return 0;
    }

	Image3dType::Pointer imageJImage = Image3dType::New();
	Image3dType::IndexType start3;
	Image3dType::SizeType  size3;
	size3[0] = xDim; size3[1] = yDim; size3[2] = zDim;
	start3[0] = 0; start3[1] = 0; start3[2] = 0;
	Image3dType::RegionType region3;
	region3.SetSize( size3 );
	region3.SetIndex( start3 );
	imageJImage->SetRegions( region3 );
	imageJImage->Allocate();
	Image3dType::PixelType value;
	value.Fill(0);
	imageJImage->FillBuffer(value);

	itk::ImageRegionIterator<Image3dType> iterRaw(imageJImage, imageJImage->GetRequestedRegion());
	// x
	for(iterRaw.GoToBegin(); !iterRaw.IsAtEnd(); ++iterRaw) {
		float x;
		rawFile >> x;
		value[0] = x;
		iterRaw.Set(value);
	}
	// y
	for(iterRaw.GoToBegin(); !iterRaw.IsAtEnd(); ++iterRaw) {
		float y;
		rawFile >> y;
		value = iterRaw.Get();
		value[1] = y;
		iterRaw.Set(value);
	}

	// translation
	Image3dType::IndexType index;
	for(iterRaw.GoToBegin(); !iterRaw.IsAtEnd(); ++iterRaw) {
		index = iterRaw.GetIndex();
		value = iterRaw.Get();
		value[0] -= index[0];
		value[1] -= index[1];
		iterRaw.Set(value);
	}

	writer->SetInput(imageJImage);
	writer->SetFileName("tumor_displacement_imagej.mhd");
	writer->Update();
	std::cout << "imagej dvf has been loaded" << std::endl;
}