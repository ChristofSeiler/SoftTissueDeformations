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

#include "Utils.h"

#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkGroupSpatialObject.h>

int main(int argc,char *argv[])
{
	// image to deform
	typedef itk::Image<unsigned char,2> GreyImageType;
	GreyImageType::Pointer image;

	// labels
	GreyImageType::Pointer label;

	// image dimension
	const int imageDim[2] = {300,300};

	// create image
	image = GreyImageType::New();
	GreyImageType::IndexType start;
	GreyImageType::SizeType  size;
	size[0] = imageDim[0]; size[1] = imageDim[1];
	start[0] = 0; start[1] = 0;
	GreyImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	image->SetRegions( region );
	image->Allocate();
	itk::ImageRegionIterator<GreyImageType> iterImage(image, image->GetRequestedRegion());

	// create label
	label = GreyImageType::New();
	label->SetRegions( region );
	label->Allocate();
	itk::ImageRegionIterator<GreyImageType> iterLabel(label, label->GetRequestedRegion());

	// save to file
	typedef itk::ImageFileWriter< GreyImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	// create image/label one tissue
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterLabel) {
		GreyImageType::IndexType index = iterImage.GetIndex();
		GreyImageType::PixelType value = 250;
		GreyImageType::PixelType label = 0;
		iterImage.Set(value);
		iterLabel.Set(label);
	}
	writer->SetInput(image);
	writer->SetFileName("OneTissue.png");
	writer->Update();
	writer->SetInput(label);
	writer->SetFileName("OneTissueLabel.png");
	writer->Update();

	// create image/label two tissues vertical
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterLabel) {
		GreyImageType::IndexType index = iterImage.GetIndex();
		GreyImageType::PixelType value;
		GreyImageType::PixelType label;
		if(index[0] >= 0 && index[0] < 150) {
			value = 150;
			label = 1;
		}
		else {
			value = 250;
			label = 0; 
		}
		iterImage.Set(value);
		iterLabel.Set(label);
	}
	writer->SetInput(image);
	writer->SetFileName("TwoTissuesVertical.png");
	writer->Update();
	writer->SetInput(label);
	writer->SetFileName("TwoTissuesVerticalLabel.png");
	writer->Update();

	// create image/label two tissues horizontal
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterLabel) {
		GreyImageType::IndexType index = iterImage.GetIndex();
		GreyImageType::PixelType value;
		GreyImageType::PixelType label;
		if(index[1] >= 0 && index[1] < 150) {
			value = 150;
			label = 1;
		}
		else {
			value = 250;
			label = 0; 
		}
		iterImage.Set(value);
		iterLabel.Set(label);
	}
	writer->SetInput(image);
	writer->SetFileName("TwoTissuesHorizontal.png");
	writer->Update();
	writer->SetInput(label);
	writer->SetFileName("TwoTissuesHorizontalLabel.png");
	writer->Update();

	// create image/label two tissues square
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterLabel) {
		GreyImageType::IndexType index = iterImage.GetIndex();
		GreyImageType::PixelType value;
		GreyImageType::PixelType label;
		if(index[1] >= 100 && index[1] < 300 && index[0] >= 0 && index[0] < 200) {
			value = 150;
			label = 1;
		}
		else {
			value = 250;
			label = 0; 
		}
		iterImage.Set(value);
		iterLabel.Set(label);
	}
	writer->SetInput(image);
	writer->SetFileName("TwoTissuesSquare.png");
	writer->Update();
	writer->SetInput(label);
	writer->SetFileName("TwoTissuesSquareLabel.png");
	writer->Update();

	// create image/label three tissues square
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterLabel) {
		GreyImageType::IndexType index = iterImage.GetIndex();
		GreyImageType::PixelType value;
		GreyImageType::PixelType label;
		if(index[1] >= 200 && index[1] < 300 && index[0] >= 0 && index[0] < 100) {
			value = 50; 
			label = 2;
		}
		else if(index[1] >= 100 && index[1] < 300 && index[0] >= 0 && index[0] < 200) {
			value = 150;
			label = 1;
		}
		else {
			value = 250; 
			label = 0; 
		}
		iterImage.Set(value);
		iterLabel.Set(label);
	}
	writer->SetInput(image);
	writer->SetFileName("ThreeTissuesSquare.png");
	writer->Update();
	writer->SetInput(label);
	writer->SetFileName("ThreeTissuesSquareLabel.png");
	writer->Update();

	// create image/label two tissues center square
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterLabel) {
		GreyImageType::IndexType index = iterImage.GetIndex();
		GreyImageType::PixelType value;
		GreyImageType::PixelType label;
		if(index[1] >= 75 && index[1] < 225 && index[0] >= 75 && index[0] < 225) {
			value = 50; 
			label = 1;
		}
		else {
			value = 150; 
			label = 0; 
		}
		iterImage.Set(value);
		iterLabel.Set(label);
	}
	writer->SetInput(image);
	writer->SetFileName("TwoTissuesCenterSquare.png");
	writer->Update();
	writer->SetInput(label);
	writer->SetFileName("TwoTissuesCenterSquareLabel.png");
	writer->Update();

	// circle with radius 75 and center at (150,150) like center square
	// create label two tissues center circle
	// create ellipse with itk
	typedef itk::EllipseSpatialObject< 2 > EllipseType;
	EllipseType::Pointer ellipse1 = EllipseType::New();
	ellipse1->SetRadius( 100.0 );
	EllipseType::TransformType::OffsetType offset;
	offset[ 0 ] = 150;
	offset[ 1 ] = 150;
	ellipse1->GetObjectToParentTransform()->SetOffset(offset);
	ellipse1->ComputeObjectToWorldTransform();

	typedef itk::SpatialObjectToImageFilter< EllipseType, GreyImageType > SpatialObjectToImageFilterType;
	SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
	imageFilter->SetInput(  ellipse1  );
	imageFilter->SetSize( size );
	SpatialObjectToImageFilterType::ValueType insideValue = 1;
	imageFilter->SetInsideValue( insideValue );
	SpatialObjectToImageFilterType::ValueType outsideValue = 0;
	imageFilter->SetOutsideValue( outsideValue );

	writer->SetInput(imageFilter->GetOutput());
	writer->SetFileName("TwoTissuesCenterCircleLabel.png");
	writer->Update();

	insideValue = 50;
	imageFilter->SetInsideValue( insideValue );
	outsideValue = 150;
	imageFilter->SetOutsideValue( outsideValue );
	writer->SetFileName("TwoTissuesCenterCircle.png");
	writer->Update();

	// circle with radius 100 and center at (150,150) like center square
	// with one tumor and one normal region
	EllipseType::Pointer ellipseSkull = EllipseType::New();
	ellipseSkull->SetRadius( 100.0 );
	offset[ 0 ] = 150;
	offset[ 1 ] = 150;
	ellipseSkull->GetObjectToParentTransform()->SetOffset(offset);
	ellipseSkull->ComputeObjectToWorldTransform();

	// tumor region
	EllipseType::Pointer ellipseTumor = EllipseType::New();
	ellipseTumor->SetRadius( 30.0 );
	offset[ 0 ] = 150;
	offset[ 1 ] = 110;
	ellipseTumor->GetObjectToParentTransform()->SetOffset(offset);
	ellipseTumor->ComputeObjectToWorldTransform();
	
	// normal region
	EllipseType::Pointer ellipseNormal = EllipseType::New();
	ellipseNormal->SetRadius( 30.0 );
	offset[ 0 ] = 150;
	offset[ 1 ] = 190;
	ellipseNormal->GetObjectToParentTransform()->SetOffset(offset);
	ellipseNormal->ComputeObjectToWorldTransform();

	typedef itk::Point<double,2> Point;
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, 
			++iterLabel) {
		GreyImageType::PixelType value = 0;
		GreyImageType::PixelType label = 250;
		Point point;
		point[0] = iterImage.GetIndex()[0];
		point[1] = iterImage.GetIndex()[1];
		if(ellipseSkull->IsInside(point)) {
			value = 1;
			label = 150;
		}
		if(ellipseTumor->IsInside(point)) {
			value = 2;
			label = 50;
		}
		if(ellipseNormal->IsInside(point)) {
			value = 0;
			label = 250;
		}
		iterImage.Set(value);
		iterLabel.Set(label);
	}
	writer->SetInput(image);
	writer->SetFileName("TumorTissuesLabel.png");
	writer->Update();
	writer->SetInput(label);
	writer->SetFileName("TumorTissues.png");
	writer->Update();

	// create grid image for deformation simulation
	const int gridSize[2] = {10,10};
	for(iterImage.GoToBegin(), iterLabel.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterLabel) {
		GreyImageType::IndexType index = iterImage.GetIndex();
		GreyImageType::PixelType value;
		int combinedIndex = index[0] * index[1];
		if(index[0] % gridSize[0] == 0 || index[1] % gridSize[1] == 0)
			value = 0;
		else
			value = 255;
		iterImage.Set(value);
	}
	writer->SetInput(image);
	writer->SetFileName("Grid.png");
	writer->Update();
}
