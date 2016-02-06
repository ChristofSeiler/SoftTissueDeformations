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

// linear interpolation of displacement in global coordinate system
void femInterpolation(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3,
					  double x4, double y4, double u1, double u2, double u3, double &u, double &v) {
	// TODO
}

int main(int argc,char *argv[])
{
	// dvf
	typedef itk::Vector<float,3> Image3dPixelType;
	typedef itk::Image<Image3dPixelType,3> Image3dType;
	typedef itk::Vector<float,2> Image2dPixelType;
	typedef itk::Image<Image2dPixelType,2> Image2dType;
	typedef itk::Image<float,3> BinaryImageType;

	// image dimension
	const int xDim = 300;
	const int yDim = 300;
	const int zDim = 1;

	// save to file
	typedef itk::ImageFileWriter< Image3dType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	typedef itk::ImageFileWriter< BinaryImageType > BinaryWriterType;
	BinaryWriterType::Pointer binaryWriter = BinaryWriterType::New();

	typedef itk::LinearInterpolateImageFunction<FloatImageType, double > InterpolatorType;
	typedef itk::InterpolateNextStepFilter<Image2dType,Image2dType,InterpolatorType> NextStepFilter;
	NextStepFilter::Pointer nextStepFilter;

	typedef itk::InterpolateNextStepFilter<Image3dType,Image3dType,InterpolatorType> Interpolateion3dFilter;
	Interpolateion3dFilter::Pointer interpolateion3dFilter;

	// load irregular fem results with itkPointSet

	typedef itk::PointSet< Image3dPixelType, 3 > PointSetType;
	PointSetType::Pointer pointSet = PointSetType::New();
	PointSetType::PointsContainer::Pointer points = PointSetType::PointsContainer::New();
	PointSetType::PointDataContainer::Pointer pointData = PointSetType::PointDataContainer::New();
	std::ifstream tumorODB("tumor_abaqus.rpy");
	if(!tumorODB) {
		std::cerr << "tumor_abaqus.rpy file not found" << std::endl;
		return 0;
	}
	std::string comment = "";
	float x[5] = {0,0,0,0,0};
	PointSetType::PointType point;
	Image3dPixelType value;
	for (unsigned int i = 0; !tumorODB.eof(); ++i) {
		tumorODB >> comment;
		tumorODB >> x[0] >> x[1] >> x[2] >> x[3] >> x[4];
		point[0] = x[0] + 150.0; point[1] = x[1] + 150.0; point[2] = 0;
		//point[0] = point[0] / 300.0 * 299.0; point[1] = point[1] / 300.0 * 299.0;
		value[0] = -x[3]; value[1] = x[4]; value[2] = 0;
		points->InsertElement( i , point );
		pointData->InsertElement( i , value );
	}
	pointSet->SetPoints( points );
	pointSet->SetPointData( pointData );

	Image3dType::Pointer tumorImage = Image3dType::New();
	Image3dType::IndexType start3;
	Image3dType::SizeType  size3;
	size3[0] = 300; size3[1] = 300; size3[2] = 1;
	start3[0] = 0; start3[1] = 0; start3[2] = 0;
	Image3dType::RegionType region3;
	region3.SetSize( size3 );
	region3.SetIndex( start3 );
	tumorImage->SetRegions( region3 );
	Image3dType::SpacingType spacing;
	//spacing[0] = (200.0)/(299.0); spacing[1] = (200.0)/(299.0); spacing[2] = 1;
	spacing[0] = 1; spacing[1] = 1; spacing[2] = 1;
	tumorImage->SetSpacing(spacing);
	tumorImage->Allocate();
	Image3dPixelType notDefined;
	notDefined.Fill(0);
	tumorImage->FillBuffer(notDefined);
	BinaryImageType::Pointer binaryImage = BinaryImageType::New();
	binaryImage->SetRegions(region3);
	binaryImage->Allocate();
	binaryImage->FillBuffer(0);

	typedef PointSetType::PointsContainer::ConstIterator  PointIterator;
	PointIterator pointItr = pointSet->GetPoints()->Begin();
	PointIterator pointEnd = pointSet->GetPoints()->End();

	typedef PointSetType::PointDataContainer::Iterator PointDataIterator;
	PointDataIterator pointDataIterator = pointData->Begin();
	
	Image3dType::IndexType index;

	while( pointItr != pointEnd ) {
		if(tumorImage->TransformPhysicalPointToIndex(pointItr.Value(),index)) {
			tumorImage->SetPixel(index, pointDataIterator.Value());
			binaryImage->SetPixel(index, 1);
		}
		pointItr++;
		pointDataIterator++;
	}
/*
	interpolateion3dFilter = Interpolateion3dFilter::New();
	interpolateion3dFilter->SetNumberOfPixels(xDim*yDim);
	interpolateion3dFilter->SetInput(tumorImage);
	interpolateion3dFilter->Update();
	tumorImage = interpolateion3dFilter->GetOutput();

	typedef itk::ResampleImageFilter<BinaryImageType,BinaryImageType> BinaryResampleFilterType;
	BinaryResampleFilterType::Pointer binaryFilter = BinaryResampleFilterType::New();
	typedef itk::AffineTransform< double, 3 > BinaryTransformType;
	BinaryTransformType::Pointer binaryTransform = BinaryTransformType::New();
	BinaryTransformType::OutputVectorType scaleb;
	scaleb[0] = (200-1.0)/(xDim-1.0);
	scaleb[1] = scaleb[0];
	scaleb[2] = 1;
	binaryTransform->Scale(scaleb);
	binaryFilter->SetTransform( binaryTransform );
	typedef itk::NearestNeighborInterpolateImageFunction<BinaryImageType, double > BinaryInterpolatorType;
	BinaryInterpolatorType::Pointer binaryInterpolator = BinaryInterpolatorType::New();
	binaryFilter->SetInterpolator( binaryInterpolator );
	double spacing3d[ 3 ];
	spacing3d[0] = 1.0; spacing3d[1] = 1.0; spacing3d[2] = 1.0;
	binaryFilter->SetOutputSpacing( spacing3d );
	double origin3d[ 3 ];
	origin3d[0] = 0.0; origin3d[1] = 0.0; origin3d[2] = 0.0;
	binaryFilter->SetOutputOrigin( origin3d );
	BinaryImageType::SizeType size3d;
	size3d[0] = tumorImage->GetBufferedRegion().GetSize()[0]; 
	size3d[1] = tumorImage->GetBufferedRegion().GetSize()[1];
	size3d[2] = 1;
	binaryFilter->SetSize(size3d);
	binaryFilter->SetInput( binaryImage );
	binaryFilter->Update();
	binaryImage = binaryFilter->GetOutput();
*/
	writer->SetInput(tumorImage);
	writer->SetFileName("tumor_displacement_abaqus_interpolated.mhd");
	writer->Update();
	binaryWriter->SetInput(binaryImage);
	binaryWriter->SetFileName("tumor_displacement_abaqus_interpolated_binary.mhd");
	binaryWriter->Update();

	std::cout << "irregular results have been loaded into point set" << std::endl;

	// load fem results with rectangular mesh

	std::ifstream rectFile("rect_abaqus.rpy");
	if(!rectFile) {
		cerr << "rect_abaqus.rpy file not found" << endl;
		return 0;
    }
	Image2dType::Pointer rectImage = Image2dType::New();
	Image2dType::IndexType start2;
	Image2dType::SizeType  size2;
	size2[0]  = 61;	size2[1]  = 61;
	start2[0] = 0; start2[1] = 0;
	Image2dType::RegionType region2;
	region2.SetSize( size2 );
	region2.SetIndex( start2 );
	rectImage->SetRegions( region2 );
	rectImage->Allocate();
	float translation[2] = {50,250};
	float shrink[2] = {0.2,0.2};
	float scale[2] = {-1,1};
	for(unsigned int i = 0; i < size2[0]*size2[1]; i++) {
		Image2dType::PixelType value;
		std::string comment;
		rectFile >> comment;
		float x[5];
		rectFile >> x[0] >> x[1] >> x[2] >> x[3] >> x[4];
		Image2dType::IndexType index;
		for(unsigned int j = 0; j < 2; j++)
			index[j] = (x[j]+translation[j])*shrink[j];
		value[0] = x[3]*scale[0];
		value[1] = x[4]*scale[1];
		rectImage->SetPixel(index, value);
	}
	nextStepFilter = NextStepFilter::New();
	nextStepFilter->SetNumberOfPixels(xDim*yDim);
	nextStepFilter->SetInput(rectImage);
	nextStepFilter->Update();

	binaryImage->FillBuffer(0);

	// convert to 3d
	Image3dType::Pointer rectParaImage = Image3dType::New();
	size3[0] = xDim; size3[1] = yDim; size3[2] = 1;
	region3.SetSize( size3 );
	region3.SetIndex( start3 );
	rectParaImage->SetRegions( region3 );
	rectParaImage->Allocate();
	itk::ImageRegionIterator<Image3dType> iterParaRect(rectParaImage, rectParaImage->GetRequestedRegion());
	itk::ImageRegionIterator<Image2dType> iterRect(nextStepFilter->GetOutput(), 
		nextStepFilter->GetOutput()->GetRequestedRegion());
	itk::ImageRegionIterator<BinaryImageType> iterBinary(binaryImage, 
		binaryImage->GetRequestedRegion());
	for(iterRect.GoToBegin(), iterParaRect.GoToBegin(), iterBinary.GoToBegin(); !iterRect.IsAtEnd(); 
			++iterRect, ++iterParaRect, ++iterBinary) {
		Image3dType::PixelType value;
		value[0] = iterRect.Get()[0];
		value[1] = iterRect.Get()[1];
		value[2] = 0.0;
		iterParaRect.Set(value);
		iterBinary.Set(1);
	}
	writer->SetInput(rectParaImage);
	writer->SetFileName("rect_displacement_abaqus_interpolated.mhd");
	writer->Update();
	binaryWriter->SetInput(binaryImage);
	binaryWriter->SetFileName("rect_displacement_abaqus_interpolated_binary.mhd");
	binaryWriter->Update();
	std::cout << "rectangular mesh has been loaded" << std::endl;
}