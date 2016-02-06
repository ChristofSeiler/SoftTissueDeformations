/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	
 */

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkInterpolateNextStepFilter.h"

int main(int argc,char *argv[])
{
	// image dimension
	const int xDim = 300;
	const int yDim = 300;
	const int zDim = 1;

	typedef itk::Vector<float,3> ImagePixelType;
	typedef itk::Image<ImagePixelType,3> ImageType;

	typedef itk::Image<float,2> FloatImageType;
	typedef itk::Image<unsigned char,2> GrayImageType;
	typedef itk::Image<float,3> BinaryImageType;

	typedef itk::ImageFileReader< ImageType > ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	typedef itk::ImageFileReader< GrayImageType > GrayReaderType;
	GrayReaderType::Pointer grayReader = GrayReaderType::New();
	typedef itk::ImageFileReader<BinaryImageType> BinaryImageReaderType;
	BinaryImageReaderType::Pointer binaryReader = BinaryImageReaderType::New();


	typedef itk::ImageFileWriter< ImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	typedef itk::ImageFileWriter< GrayImageType > GrayWriterType;
	GrayWriterType::Pointer grayWriter = GrayWriterType::New();
	typedef itk::ImageFileWriter< BinaryImageType > BinaryWriterType;
	BinaryWriterType::Pointer binaryWriter = BinaryWriterType::New();

	typedef itk::NearestNeighborInterpolateImageFunction<FloatImageType, double > InterpolatorType;
	typedef itk::InterpolateNextStepFilter<ImageType,ImageType,InterpolatorType > ImageFilterType;
	ImageFilterType::Pointer scaleFilter;

	reader->SetFileName("input/TwoDirectionExpansion.mhd");

	scaleFilter = ImageFilterType::New();
	scaleFilter->SetNumberOfPixels(xDim*yDim);
	scaleFilter->SetInput(reader->GetOutput());

	writer->SetInput(scaleFilter->GetOutput());
	writer->SetFileName("output/TwoDirectionExpansion.mhd");
	writer->Update();

	reader->SetFileName("input/TwoDirectionExpansionConfidence.mhd");
	scaleFilter->SetInput(reader->GetOutput());
	writer->SetInput(scaleFilter->GetOutput());
	writer->SetFileName("output/TwoDirectionExpansionConfidence.mhd");
	writer->Update();

	reader->SetFileName("input/rect_displacement_abaqus_interpolated.mhd");
	scaleFilter->SetInput(reader->GetOutput());
	writer->SetInput(scaleFilter->GetOutput());
	writer->SetFileName("output/rect_displacement_abaqus_interpolated.mhd");
	writer->Update();

	binaryReader->SetFileName("input/rect_displacement_abaqus_interpolated_binary.mhd");
	binaryReader->Update();

	typedef itk::ResampleImageFilter<BinaryImageType,BinaryImageType> BinaryResampleFilterType;
	BinaryResampleFilterType::Pointer binaryFilter = BinaryResampleFilterType::New();
	typedef itk::AffineTransform< double, 3 > BinaryTransformType;
	BinaryTransformType::Pointer binaryTransform = BinaryTransformType::New();
	BinaryTransformType::OutputVectorType scale;
	scale[0] = (binaryReader->GetOutput()->GetBufferedRegion().GetSize()[0]-1.0)/(xDim-1.0);
	scale[1] = scale[0];
	scale[2] = 1;
	binaryTransform->Scale(scale);
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
	size3d[0] = xDim; size3d[1] = yDim; size3d[2] = 1;
	binaryFilter->SetSize(size3d);
	binaryFilter->SetInput( binaryReader->GetOutput() );

	binaryWriter->SetInput(binaryFilter->GetOutput());
	binaryWriter->SetFileName("output/rect_displacement_abaqus_interpolated_binary.mhd");
	binaryWriter->Update();

	grayReader->SetFileName("input/label_image.png");
	grayReader->Update();

	typedef itk::ResampleImageFilter<GrayImageType,GrayImageType> ResampleFilterType;
	ResampleFilterType::Pointer grayFilter = ResampleFilterType::New();
	typedef itk::AffineTransform< double, 2 > TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->Scale((grayReader->GetOutput()->GetBufferedRegion().GetSize()[0]-1.0)/(xDim-1.0));
	grayFilter->SetTransform( transform );
	typedef itk::NearestNeighborInterpolateImageFunction<GrayImageType, double > GrayInterpolatorType;
	GrayInterpolatorType::Pointer interpolator = GrayInterpolatorType::New();
	grayFilter->SetInterpolator( interpolator );
	double spacing[ 2 ];
	spacing[0] = 1.0; spacing[1] = 1.0;
	grayFilter->SetOutputSpacing( spacing );
	double origin[ 2 ];
	origin[0] = 0.0; origin[1] = 0.0;
	grayFilter->SetOutputOrigin( origin );
	GrayImageType::SizeType size;
	size.Fill(xDim);
	grayFilter->SetSize(size);
	grayFilter->SetInput( grayReader->GetOutput() );

	grayWriter->SetInput(grayFilter->GetOutput());
	grayWriter->SetFileName("output/label_image.png");
	grayWriter->Update();

	grayReader->SetFileName("input/tissue_image.png");
	grayFilter->SetInput(grayReader->GetOutput());
	grayWriter->SetInput(grayFilter->GetOutput());
	grayWriter->SetFileName("output/tissue_image.png");
	grayWriter->Update();

	grayReader->SetFileName("input/Grid.png");
	grayFilter->SetInput(grayReader->GetOutput());
	grayWriter->SetInput(grayFilter->GetOutput());
	grayWriter->SetFileName("output/Grid.png");
	grayWriter->Update();
}
