/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	For general definitions.
 */

#pragma once

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkImageToVTKImageFilter.h>
#include <itkVTKImageToImageFilter.h>

// type definitions itk
typedef float PixelType;
const unsigned int Dimension = 2;
typedef itk::Image<PixelType,Dimension> ImageType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::ImageToVTKImageFilter<ImageType> VTKExporterType;
typedef itk::VTKImageToImageFilter<ImageType> ITKExporterType;
typedef itk::ImageFileWriter<ImageType> ImageWriterType;
typedef itk::Image<unsigned char,2> GreyImageType;
const unsigned int TissueTypesMax = 10;
typedef itk::Vector<float,TissueTypesMax> NLabelPixelType;
typedef itk::Image<NLabelPixelType,2> NLabelImageType;

// vector
typedef itk::Vector<float,Dimension> VectorPixelType;
typedef itk::Image<VectorPixelType,Dimension> VectorImageType;
typedef itk::ImageToVTKImageFilter<VectorImageType> VectorVTKExporterType;
typedef itk::VTKImageToImageFilter<VectorImageType> VectorITKExporterType;
typedef itk::ImageFileWriter<VectorImageType> VectorImageWriterType;

typedef itk::Vector<float,3> Vector3PixelType;
typedef itk::Image<Vector3PixelType,3> Vector3ImageType;

typedef itk::Image<float,2> FloatImageType;
typedef itk::Image<float,3> BinaryImageType;
