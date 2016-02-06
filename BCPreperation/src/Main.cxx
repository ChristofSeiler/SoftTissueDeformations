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
#include <itkMRFRegularizationFilterGPU.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileReader.h>
#include <itkEllipseSpatialObject.h>

int main(int argc,char *argv[])
{
	// dvf
	typedef itk::Vector<float,3> ImagePixelType;
	typedef itk::Image<ImagePixelType,3> ImageType;

	// image dimension
	const int xDim = 300;
	const int yDim = 300;
	const int zDim = 1;

	// create image
	ImageType::Pointer image = ImageType::New();
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0] = xDim; size[1] = yDim; size[2] = zDim;
	start[0] = 0; start[1] = 0; start[2] = 0;
	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );
	image->SetRegions( region );
	image->Allocate();

	// create float image
	ImageType::Pointer confidenceImage = ImageType::New();
	confidenceImage->SetRegions(region);
	confidenceImage->Allocate();

	// save to file
	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();

	// sphincter like expansion

	const float r = 100.0;
	const float a = 150.0;
	const float b = 150.0;
	const float magnitude = 50.0;

	itk::ImageRegionIterator<ImageType> iter(image, image->GetRequestedRegion());
	itk::ImageRegionIterator<ImageType> confidenceIter(confidenceImage, confidenceImage->GetRequestedRegion());

	for(iter.GoToBegin(), confidenceIter.GoToBegin(); !iter.IsAtEnd(); ++iter, ++confidenceIter) {
		ImageType::PixelType value;
		value[0] = 0; value[1] = 0; value[2] = 0;
		float halfx = image->GetRequestedRegion().GetSize()[0]/2;
		float halfy = image->GetRequestedRegion().GetSize()[1]/2;
		if(iter.GetIndex()[0] < halfx) {
			value[0] = (halfx-iter.GetIndex()[0])/halfx*magnitude;
		}
		else {
			value[0] = -(iter.GetIndex()[0]+1-halfx)/halfx*magnitude;
		}
		if(iter.GetIndex()[1] < halfy) {
			value[1] = (halfy-iter.GetIndex()[1])/halfy*magnitude;
		}
		else {
			value[1] = -(iter.GetIndex()[1]+1-halfy)/halfy*magnitude;
		}
		iter.Set(value);
	}

	// create ellipse with itk
	typedef itk::EllipseSpatialObject< 3 > EllipseType;
	EllipseType::Pointer ellipse1 = EllipseType::New();
	ellipse1->SetRadius( r );
	EllipseType::TransformType::OffsetType offset;
	offset[ 0 ] = a;
	offset[ 1 ] = b;
	offset[ 2 ] = 0.0;
	ellipse1->GetObjectToParentTransform()->SetOffset(offset);
	ellipse1->ComputeObjectToWorldTransform();

	EllipseType::PointType point;
	for(confidenceIter.GoToBegin(); !confidenceIter.IsAtEnd(); ++confidenceIter) {
		ImageType::PixelType confidenceValue;
		confidenceValue.Fill(1);
		point[0] = confidenceIter.GetIndex()[0];
		point[1] = confidenceIter.GetIndex()[1]; 
		point[2] = confidenceIter.GetIndex()[2];
		if(ellipse1->IsInside(point))
			confidenceValue.Fill(0);
		confidenceIter.Set(confidenceValue);
	}

	writer->SetInput(image);
	writer->SetFileName("SphincterExpansion.mhd");
	writer->Update();
	writer->SetInput(confidenceImage);
	writer->SetFileName("SphincterExpansionConfidence.mhd");
	writer->Update();

	// sqaure expansion

	for(iter.GoToBegin(), confidenceIter.GoToBegin(); !iter.IsAtEnd(); ++iter, ++confidenceIter) {
		ImageType::PixelType value;
		ImageType::PixelType confidenceValue;
		confidenceValue[0] = 0; confidenceValue[1] = 0; confidenceValue[2] = 0;
		value[0] = 0; value[1] = 0; value[2] = 0;
		if(iter.GetIndex()[0] >= 299 && iter.GetIndex()[0] <= 299) {
			value[0] = -magnitude;
			confidenceValue[0] = 1;
		}
		else if(iter.GetIndex()[0] >= 0 && iter.GetIndex()[0] <= 0) {
			value[0] = magnitude;
			confidenceValue[0] = 1;
		}
		if(iter.GetIndex()[1] >= 299 && iter.GetIndex()[1] <= 299) {
			value[1] = -magnitude;
			confidenceValue[1] = 1;
		}
		else if(iter.GetIndex()[1] >= 0 && iter.GetIndex()[1] <= 0) {
			value[1] = magnitude;
			confidenceValue[1] = 1;
		}
		iter.Set(value);
		confidenceIter.Set(confidenceValue);
	}
	writer->SetInput(image);
	writer->SetFileName("SquareExpansion.mhd");
	writer->Update();
	writer->SetInput(confidenceImage);
	writer->SetFileName("SquareExpansionConfidence.mhd");
	writer->Update();

	// two direction expansion

	for(iter.GoToBegin(), confidenceIter.GoToBegin(); !iter.IsAtEnd(); ++iter, ++confidenceIter) {
		ImageType::PixelType value;
		ImageType::PixelType confidenceValue;
		confidenceValue[0] = 0; confidenceValue[1] = 0; confidenceValue[2] = 0;
		value[0] = 0; value[1] = 0; value[2] = 0;
		if(iter.GetIndex()[0] >= 299 && iter.GetIndex()[0] <= 299) {
			value[0] = -50;
			confidenceValue[0] = 1;
		}
		else if(iter.GetIndex()[0] >= 0 && iter.GetIndex()[0] <= 0) {
			value[0] = 0;
			confidenceValue[0] = 1;
		}
		if(iter.GetIndex()[1] >= 0 && iter.GetIndex()[1] <= 0) {
			value[1] = 50;
			confidenceValue[1] = 1;
		}
		else if(iter.GetIndex()[1] >= 299 && iter.GetIndex()[1] <= 299) {
			value[1] = 0;
			confidenceValue[1] = 1;
		}
		iter.Set(value);
		confidenceIter.Set(confidenceValue);
	}
	writer->SetInput(image);
	writer->SetFileName("TwoDirectionExpansion.mhd");
	writer->Update();
	writer->SetInput(confidenceImage);
	writer->SetFileName("TwoDirectionExpansionConfidence.mhd");
	writer->Update();

}
