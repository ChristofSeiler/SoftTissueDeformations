/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Utils for general use.
 */

#include "Utils.h"

#include "MechanicalProperties.h"

#include <vtkImageData.h>
#include <vtkImageViewer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageActor.h>
#include <vtkImageViewer.h>
#include <vtkActor2D.h>
#include <vtkGridTransform.h>
#include <vtkImageReslice.h>
#include <vtkCellPicker.h>
#include <vtkImagePlaneWidget.h>
#include <vtkProperty.h>
#include <vtkTextProperty.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>
#include <vtkWindowLevelLookupTable.h>
#include <vtkGlyphSource2D.h>
#include <vtkGlyph2D.h>
#include <vtkArrowSource.h>
#include <vtkMaskPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkExtractGrid.h>
#include <vtkImageShrink3D.h>
#include <vtkPolyData.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>

#include <itkVector.h>
#include <itkImage.h>
#include <itkImageToVTKImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkResampleImageFilter.h>
#include <itkNeighborhoodIterator.h>
#include <itkConstantBoundaryCondition.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include "DisplacementBoundaryCondition.h"

using namespace dvfRegularization;
using namespace std;

void Utils::visualizeImage(vtkImageData *im) {
	vtkImageViewer * viewer = vtkImageViewer::New();
	vtkRenderWindowInteractor * renderWindowInteractor = vtkRenderWindowInteractor::New();
	viewer->SetupInteractor( renderWindowInteractor );
	viewer->SetInput( im );
	viewer->SetColorWindow( 256 );
	viewer->SetColorLevel( 128 );
	viewer->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
	viewer->Render();
	renderWindowInteractor->Start();
	viewer->Delete();
	renderWindowInteractor->Delete();
}

void Utils::visualizeMultipleImages(vector<vtkImageData*> &images) {
	// create component
	vtkRenderWindowInteractor *renderWindowInteractor = vtkRenderWindowInteractor::New();
	vtkRenderWindow *renderWindow = vtkRenderWindow::New();

	int width = 0;
	int height = 0;
	for(unsigned int i = 0; i < images.size(); i++) {
		vtkRenderer *renderer = vtkRenderer::New();
		vtkImageMapper *imageMapper  = vtkImageMapper::New();
		vtkActor2D *actor2D = vtkActor2D::New();

		// setup the pipeline
		actor2D->SetMapper(imageMapper);
		renderer->AddActor2D(actor2D);
		renderWindow->AddRenderer(renderer);

		// set parameter
		imageMapper->SetInput(images[i]);
		imageMapper->SetColorWindow(200);
		imageMapper->SetColorLevel(0);
		renderer->SetViewport(i*(1.0/images.size()), 0, 1.0/images.size()+i*(1.0/images.size()), 1);
		renderer->SetBackground(1.0, 1.0, 1.0);

		int *ext = imageMapper->GetInput()->GetWholeExtent();
		width = ext[1] - ext[0] + 1;
		height = ext[3] - ext[2] + 1;
	}

	// window size
	renderWindow->SetSize(images.size()*width+50, height+50);

	// render
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindow->Render();
	renderWindowInteractor->Start();

	// delete
	renderWindow->Delete();
	renderWindowInteractor->Delete();
}

void Utils::visualizeMultipleImages3D(vector<vtkImageData*> &images) {
	vtkRenderWindow *renderWin=vtkRenderWindow::New();
	int *screen = renderWin->GetScreenSize();
	int width = *screen;
	int height = *(screen+1);
	renderWin->SetSize(3*(0.5*height), 0.5*height);
	for(unsigned int i = 0; i < images.size(); i++) {
		// vtkImageActor is used to render an image in a 3D scene
		vtkImageActor *actor=vtkImageActor::New();
		actor->SetInput(images[i]);
  		int ext[6];
 		actor->GetInput()->GetExtent( ext );
		ext[4]=0;ext[5]=0;
		actor->SetDisplayExtent( ext );

		vtkRenderer *renderer=vtkRenderer::New();
		renderer->SetBackground(1.0, 1.0, 1.0);
		renderer->SetViewport(i*(1.0/images.size()), 0, 1.0/images.size()+i*(1.0/images.size()), 1);
		renderer->AddActor(actor);
		
		renderWin->AddRenderer(renderer);
	}

	// Create an interactor window (play with actor (object) by mouse)
	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	// Add the render window to the interactor
	iren->SetRenderWindow(renderWin);

	renderWin->Render();
	iren->Start();

	// delete
	renderWin->Delete();
	iren->Delete();
}

void Utils::visualizeImage3D(vtkImageData *im) {
	vtkRenderer *renderer=vtkRenderer::New();
	renderer->SetBackground(1,1,1);

	vtkRenderWindow* renWin = vtkRenderWindow::New();
	vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();

	renWin->AddRenderer(renderer);
	iren->SetRenderWindow(renWin);

	vtkCellPicker* picker = vtkCellPicker::New();
	picker->SetTolerance(0.005);
	vtkImagePlaneWidget *iwidget=vtkImagePlaneWidget::New();
	iwidget->SetInteractor(iren);
	iwidget->SetKeyPressActivationValue('x');
	iwidget->SetPicker(picker);
	iwidget->RestrictPlaneToVolumeOn();
	iwidget->GetPlaneProperty()->SetColor(1,0,0);
	iwidget->GetPlaneProperty()->SetOpacity(0.5);
	iwidget->GetTextProperty()->SetColor(0,0,1);
	iwidget->SetInput(im);
	iwidget->SetPlaneOrientationToZAxes();
	iwidget->SetSliceIndex(10);
	iwidget->DisplayTextOn();
	//iwidget->GetImageMapToColors()->SetOuputFormatToRGB();
	// iwidget->GetImageMapToColors()->GetLookupTable()->SetAlpha(0.8);
	// iwidget->On();
	//planeWidgetX->InteractionOff();
	//planeWidgetX->InteractionOn();

	renWin->Render();
	iren->Start();

	renderer->Delete();
	picker->Delete();
	iwidget->Delete();
	renWin->Delete();
	iren->Delete();
}

void Utils::visualizeDvfArrow(vtkPolyData* image, const float factor) {
	vtkArrowSource *arrow = vtkArrowSource::New();
	vtkGlyph3D *glyph = vtkGlyph3D::New();
	glyph->SetInput( image );
	glyph->SetSourceConnection( arrow->GetOutputPort() );
	glyph->SetScaleModeToScaleByScalar();
	glyph->SetColorModeToColorByScalar();
	glyph->SetVectorModeToUseVector();
	glyph->OrientOn();
	glyph->SetScaleFactor(factor);
	vtkPolyDataMapper *glyphmapper = vtkPolyDataMapper::New();
	glyphmapper->SetInputConnection(glyph->GetOutputPort());
	vtkActor *glyphactor = vtkActor::New();
	glyphactor->SetMapper(glyphmapper);
	vtkRenderer *ren1 = vtkRenderer::New();
	ren1->AddActor(glyphactor);

	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
	vtkRenderWindow *renWin = vtkRenderWindow::New();
	renWin->AddRenderer(ren1);
	iren->SetRenderWindow(renWin);
	renWin->Render();
	iren->Start();
}

void Utils::compoundImage(vtkImageData *dvfx,vtkImageData *dvfy,vtkImageData *dvfz,vtkImageData *DVF)
{
	typedef float dvf_scalar;
	dvf_scalar *px,*py,*pz,*Px,*Py,*Pz;
	int dim[3];
	dvfx->GetDimensions(dim);

	for ( int i = 0 ; i < dim[0] ; i++)
		for ( int j = 0 ; j < dim[1] ; j++)
			for ( int k = 0 ; k < dim[2] ; k++)
			{
				px=(dvf_scalar*)dvfx->GetScalarPointer(i,j,k);
				py=(dvf_scalar*)dvfy->GetScalarPointer(i,j,k);
				pz=(dvf_scalar*)dvfz->GetScalarPointer(i,j,k);
			
				Px=(dvf_scalar*)DVF->GetScalarPointer(i,j,k);

				*Px=*px;
				Py=Px+1;
				Pz=Py+1;

				*Py=*py;
				*Pz=*pz;
			}
}

void Utils::deformImage(vtkImageData *inputImage,vtkImageData *DVF,vtkImageData *outputImage)
{
	vtkGridTransform *gridTransform=vtkGridTransform::New();
	gridTransform->SetDisplacementGrid(DVF);
	//gridTransform->SetDisplacementShift(5);
	//gridTransform->SetDisplacementShift(grid->GetDisplacementShift());
	//gridTransform->SetDisplacementScale(grid->GetDisplacementScale());

	gridTransform->Update();

	vtkImageReslice *reslice=vtkImageReslice::New();
	reslice->SetInput(inputImage);
	reslice->SetResliceTransform(gridTransform);
	reslice->SetInterpolationModeToNearestNeighbor();
	reslice->SetOutputSpacing(inputImage->GetSpacing());

	reslice->Update();

	outputImage->DeepCopy(reslice->GetOutput());

	gridTransform->Delete();
	reslice->Delete();
}

vtkImageData* Utils::scalarToColor(vtkImageData* image) {
	// color mapping
	// Create a greyscale lookup table
	vtkWindowLevelLookupTable *lut = vtkWindowLevelLookupTable::New();
	lut->SetNumberOfTableValues(2000);
	int index = 0;
	/*
	for(int i = lut->GetNumberOfTableValues()/2; i > 0; i--) {
		lut->SetTableValue( index, i*1.0/(lut->GetNumberOfTableValues()/2), 0.0, 0.0, 1.0 );
		index++;
	}
	for(int i = 0; i < lut->GetNumberOfTableValues()/2; i++) {
		lut->SetTableValue( index, 0.0, 0.0, i*1.0/(lut->GetNumberOfTableValues()/2), 1.0 );
		index++;
	}*/
	for(int i = 0; i < lut->GetNumberOfTableValues(); i++) {
		lut->SetTableValue( index, 0.0, 0.0, i*1.0/(lut->GetNumberOfTableValues()/2), 1.0 );
		index++;
	}
	lut->SetLevel(0);
	lut->SetWindow(100);
	lut->Build();

	// Map the image through the lookup table
	vtkImageMapToColors *color = vtkImageMapToColors::New();
	color->SetLookupTable(lut);

	color->SetInput(image);
	color->Update();
	return color->GetOutput();
}

vtkPolyData* Utils::convertImageToPoly(vtkImageData* image) {
	// Create the dataset. In this case, we create a vtkPolyData
	vtkPolyData* polydata = vtkPolyData::New();
	// Assign points
	vtkPoints *points = vtkPoints::New();

	vtkFloatArray *scalars = vtkFloatArray::New();
	scalars->SetNumberOfComponents(1);
	vtkFloatArray *vectors = vtkFloatArray::New();
	vectors->SetNumberOfComponents(3);

	int dims[3];
	image->GetDimensions(dims);
	int subset[3] = { 10, 10, 1 };
	int noOfTubles = dims[0]/subset[0]*dims[1]/subset[1]*dims[2]/subset[2];
	scalars->SetNumberOfTuples(noOfTubles);
	vectors->SetNumberOfTuples(noOfTubles);
	points->Allocate(noOfTubles);

	double spacing[3];
	image->GetSpacing(spacing);
	int id = 0;

	for ( int k=0; k<dims[2]; k += subset[2]) {
		for ( int j=0; j<dims[1]; j += subset[1]) {
			for ( int i=0; i<dims[0]; i += subset[0]) {
				float *currentPixel = (float*)image->GetScalarPointer(i, j, k);
				float v[3];
				v[0] = *currentPixel; currentPixel++;
				v[1] = *currentPixel; currentPixel++;
				v[2] = *currentPixel;
				points->InsertPoint(id, i*spacing[0], j*spacing[1], k*spacing[2]);
				float magnitude[1] = { sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) };
				scalars->InsertTuple(id, magnitude);
				vectors->InsertTuple(id, v);
				id++;
			}
		}
	}
	polydata->SetPoints(points);
	// Assign scalars
	polydata->GetPointData()->SetScalars(scalars);
	// Assign vectors
	polydata->GetPointData()->SetVectors(vectors);
	return polydata;
}

vtkPolyData* Utils::convertImageToPoly(const VectorImageType::Pointer image) {
	// Create the dataset. In this case, we create a vtkPolyData
	vtkPolyData* polydata = vtkPolyData::New();
	// Assign points
	vtkPoints *points = vtkPoints::New();

	vtkFloatArray *scalars = vtkFloatArray::New();
	scalars->SetNumberOfComponents(1);
	vtkFloatArray *vectors = vtkFloatArray::New();
	vectors->SetNumberOfComponents(3);

	VectorImageType::RegionType region = image->GetRequestedRegion();
	int noOfTubles = region.GetSize()[0]*region.GetSize()[1]*region.GetSize()[2];
	scalars->SetNumberOfTuples(noOfTubles);
	vectors->SetNumberOfTuples(noOfTubles);
	points->Allocate(noOfTubles);

	double spacing[3] = {1,1,1};
	int id = 0;

	itk::ImageRegionConstIterator<VectorImageType> it(image, image->GetRequestedRegion());
	for(; !it.IsAtEnd(); ++it) {
		float v[3];
		v[0] = it.Get()[0];
		v[1] = it.Get()[1];
		v[2] = it.Get()[2];
		points->InsertPoint(id, it.GetIndex()[0]*spacing[0],
			it.GetIndex()[1]*spacing[1], 
			it.GetIndex()[2]*spacing[2]);
		float magnitude[1] = { sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) };
		scalars->InsertTuple(id, magnitude);
		vectors->InsertTuple(id, v);
		id++;		
	}
	polydata->SetPoints(points);
	// Assign scalars
	polydata->GetPointData()->SetScalars(scalars);
	// Assign vectors
	polydata->GetPointData()->SetVectors(vectors);
	return polydata;
}

VectorImageType::Pointer Utils::convertPolyToITK(vtkPolyData* poly) {
	VectorImageType::Pointer image = VectorImageType::New();

	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = 21;	size[1]  = 21; size[2] = 1;
	start[0] = 0; start[1] = 0; start[2] = 0;

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	// iterate
	itk::ImageRegionIterator<VectorImageType> iter(image, image->GetRequestedRegion());
	vtkFloatArray* vectors = dynamic_cast<vtkFloatArray*>(poly->GetPointData()->GetVectors());
	for(unsigned int i = 0; !iter.IsAtEnd(); ++iter, ++i) {
		VectorImageType::PixelType pixelValue;
		double* v = vectors[i].GetTuple(i);
		pixelValue[0] = v[0];
		pixelValue[1] = v[1];
		pixelValue[2] = v[2];
		iter.Set(pixelValue);
	}

	return image;
}

VectorImageType::Pointer Utils::createVectorImage(ImageType::Pointer elastic[3]) {
	VectorImageType::Pointer image = VectorImageType::New();
	ImageType::SizeType oneSize = elastic[0]->GetRequestedRegion().GetSize();
	ImageType::IndexType oneIndex = elastic[0]->GetRequestedRegion().GetIndex();

	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = oneSize[0];  // size along X
	size[1]  = oneSize[1];  // size along Y
	size[2] = 1;

	start[0] = oneIndex[0];  // first index on X
	start[1] = oneIndex[1];  // first index on Y
	start[2] = 0;

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	// iterate
	itk::ImageRegionIterator<ImageType> e0(elastic[0], elastic[0]->GetRequestedRegion());
	itk::ImageRegionIterator<ImageType> e1(elastic[1], elastic[1]->GetRequestedRegion());
	itk::ImageRegionIterator<ImageType> e2(elastic[2], elastic[2]->GetRequestedRegion());
	itk::ImageRegionIterator<VectorImageType> iterImage(image, image->GetRequestedRegion());
	for(; !iterImage.IsAtEnd(); ++iterImage, ++e0, ++e1, ++e2) {
		VectorImageType::PixelType pixelValue;
		pixelValue[0] = e0.Get();
		pixelValue[1] = e1.Get();
		pixelValue[2] = e2.Get();
		iterImage.Set(pixelValue);
	}

	return image;
}

vtkImageData* Utils::convertSampledVectorImage(VectorImageType::Pointer image) {
	vtkImageData* convertedImage = vtkImageData::New();
	// input image vtkImageData's
	int dim[3] = {21,21,1};
	//inputImage->SetExtent( 0, width, 0, height, 0, 0);
	convertedImage->SetDimensions(dim);
	convertedImage->SetOrigin( 0.0, 0.0, 0.0);
	convertedImage->SetSpacing( 1.0, 1.0, 1.0);
	convertedImage->SetScalarTypeToFloat();
	convertedImage->SetNumberOfScalarComponents(3);
	convertedImage->AllocateScalars();
	itk::ImageRegionIterator<VectorImageType> iterImage(image, image->GetRequestedRegion());
	for(int j = 0; j < dim[1]; j++)
	{
		for(int i = 0; i < dim[0]; i++)
		{
			float *currentPixel = (float*) convertedImage->GetScalarPointer( i , j , 0 );
			VectorImageType::PixelType pixelValue = iterImage.Value();
			*currentPixel = pixelValue[0]; currentPixel++;
			*currentPixel = pixelValue[1]; currentPixel++;
			*currentPixel = pixelValue[2];
			++iterImage;
		}
	}
	return convertedImage;
}

vtkImageData* Utils::convertVectorImage(VectorImageType::Pointer image) {
	vtkImageData* convertedImage = vtkImageData::New();
	// input image vtkImageData's
	VectorImageType::SizeType size = image->GetRequestedRegion().GetSize();
	int dim[3] = {size[0],size[1],size[2]};
	//inputImage->SetExtent( 0, width, 0, height, 0, 0);
	convertedImage->SetDimensions(dim);
	convertedImage->SetOrigin( 0.0, 0.0, 0.0);
	convertedImage->SetSpacing( 1.0, 1.0, 1.0);
	convertedImage->SetScalarTypeToFloat();
	convertedImage->SetNumberOfScalarComponents(3);
	convertedImage->AllocateScalars();
	itk::ImageRegionIterator<VectorImageType> iterImage(image, image->GetRequestedRegion());
	for(int j = 0; j < dim[1]; j++)
	{
		for(int i = 0; i < dim[0]; i++)
		{
			float *currentPixel = (float*) convertedImage->GetScalarPointer( i , j , 0 );
			VectorImageType::PixelType pixelValue = iterImage.Value();
			*currentPixel = pixelValue[0]; currentPixel++;
			*currentPixel = pixelValue[1]; currentPixel++;
			*currentPixel = pixelValue[2];
			++iterImage;
		}
	}
	return convertedImage;
}

void Utils::createGridLineImage(vtkImageData* image, int gridSize[]) {
	image->SetOrigin( 0.0, 0.0, 0.0);
	image->SetSpacing( 1.0, 1.0, 1.0);
	image->SetScalarTypeToFloat();
	image->SetNumberOfScalarComponents(1);
	image->AllocateScalars();
	int dim[3];
	image->GetDimensions(dim);
	for(int i = 0; i < dim[0]; i++)
	{
		for(int j = 0; j < dim[1]; j++)
		{
			int index = i * j;
			// get scalar pointer to current pixel
			if(i % gridSize[0] == 0 || j % gridSize[1] == 0) {
				float *currentPixel = (float*) image->GetScalarPointer( i , j , 0 );
				*currentPixel = 80;
			}
		}
	}
}

VectorImageType::Pointer Utils::createComplicatedInitialGuessDVF(const int imageSize[], 
										const float maxDeform, 
										const ImageType::Pointer labels, 
										const std::vector<MechanicalProperties*> &mechanicalProperties) {
	VectorImageType::Pointer image = VectorImageType::New();
	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	float total = 0;
	for(int i = 0; i < mechanicalProperties.size(); i++)
		total += mechanicalProperties[i]->getYoungsModulus();

	itk::NeighborhoodIterator<VectorImageType>::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;
	typedef itk::ConstantBoundaryCondition<VectorImageType> BC;
	itk::NeighborhoodIterator<VectorImageType,BC> imageIter(radius, image, image->GetRequestedRegion());

	// count tisse pixel per vertical and horizontal line
	itk::ImageRegionConstIterator<ImageType> labelsIter(labels, labels->GetRequestedRegion());
	std::vector<Count*> allTissues(mechanicalProperties.size());
	// initialize vectors
	for(int i = 0; i < allTissues.size(); i++)
		allTissues[i] = new Count(new std::vector<int>(imageSize[1]), new std::vector<int>(imageSize[0]));
	// counting
	for(imageIter.GoToBegin(); !labelsIter.IsAtEnd(); ++labelsIter) {
		ImageType::IndexType index = labelsIter.GetIndex();
		ImageType::PixelType pixel = labelsIter.Get();
		allTissues[pixel]->horizontalLine->at(index[1])++;
		allTissues[pixel]->verticalLine->at(index[0])++;
	}

	/*
	std::cout << "horizontal line:" << std::endl;
	for(int i = 0; i < allTissues[0]->horizontalLine->size(); i++)
		std::cout << "position " << i << " t0 " << allTissues[0]->horizontalLine->at(i) << " t1 " << 
			allTissues[1]->horizontalLine->at(i) << " t2 " << allTissues[2]->horizontalLine->at(i) << std::endl;

	std::cout << "vertical line:" << std::endl;
	for(int i = 0; i < allTissues[0]->verticalLine->size(); i++)
		std::cout << "position " << i << " t0 " << allTissues[0]->verticalLine->at(i) << " t1 " << 
			allTissues[1]->verticalLine->at(i) << " t2 " << allTissues[2]->verticalLine->at(i) << std::endl;
	*/

	for(imageIter.GoToBegin(), labelsIter.GoToBegin(); !imageIter.IsAtEnd(); ++imageIter, ++labelsIter) {
		// get pixel property
		ImageType::PixelType type = labelsIter.Get();
		ImageType::PixelType ym = mechanicalProperties[type]->getYoungsModulus();
		ImageType::PixelType pr = mechanicalProperties[type]->getPoissonRatio();
		// index of pixel
		VectorImageType::IndexType index = imageIter.GetIndex();

		// calc ideal change in displacement in x direction
		ImageType::PixelType ratio = ym/total;
		ImageType::PixelType globalChangeX = (1-ratio)*maxDeform;
		ImageType::PixelType localChangeX = globalChangeX/allTissues[type]->horizontalLine->at(index[1]);

		// calc ideal change in displacement in y direction
		//float nonAir = imageSize[0]-allTissues[0]->horizontalLine->at(index[1]);
		//float nonAirRatio = nonAir/imageSize[0];
		ImageType::PixelType globalChangeY = globalChangeX*pr/2;
		ImageType::PixelType localChangeY = globalChangeY/(allTissues[type]->verticalLine->at(index[0])/2);

		// get displacement of previous pixel
		ImageType::PixelType prevDispX = imageIter.GetPixel(3)[0];
		ImageType::PixelType prevDispY = imageIter.GetPixel(1)[1];

		// change sign in x
		int sign = -1;

		// continue deformation
		VectorImageType::PixelType elastic;
		elastic[0] = prevDispX+sign*localChangeX;
		if(index[0] == 200)
			std::cout << "elastic[0] " << elastic[0] << std::endl;

		// boundary condition
		if(index[1] == 0)
			elastic[1] = sign*globalChangeY;
		else
			elastic[1] = prevDispY+localChangeY;

		elastic[2] = 0;
		imageIter.SetCenterPixel(elastic);
	}

	// delete resources
	for(int i = 0; i < allTissues.size(); i++) delete allTissues[i];

	return image;
}

VectorImageType::Pointer Utils::createEmtpyDVF(int imageSize[]) {
	VectorImageType::Pointer image = VectorImageType::New();
	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();
	
	itk::ImageRegionIterator<VectorImageType> iter(image, image->GetRequestedRegion());
	for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
		VectorImageType::PixelType value;
		value[0] = 0; value[1] = 0; value[2] = 0;
		iter.Set(value);
	}

	return image;
}

VectorImageType::Pointer Utils::extrapolateNearestNeighborInitialDVF(int pixelCount, 
																	 VectorImageType::Pointer field) {
	VectorImageType::Pointer image = VectorImageType::New();
	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = std::sqrt((float)pixelCount);  // size along X
	size[1]  = std::sqrt((float)pixelCount);  // size along Y
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	if(pixelCount == 1) {
		VectorImageType::IndexType index;
		index[0] = 0; index[1] = 0;
		VectorImageType::PixelType value;
		value[0] = 0; value[1] = 0;
		image->SetPixel(index, value);
	}
	else {
		itk::ImageRegionIterator<VectorImageType> iter(image, image->GetRequestedRegion());
		for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
			VectorImageType::IndexType index = iter.GetIndex();
			float smaller = field->GetRequestedRegion().GetSize()[0];
			float bigger = image->GetRequestedRegion().GetSize()[0];
			float ratio = smaller / bigger;
			index[0] = index[0] * ratio;
			index[1] = index[1] * ratio;
			iter.Set(field->GetPixel(index));
		}
	}

	return image;
}

VectorImageType::Pointer Utils::extrapolateLinearInitialDVF(int pixelCount, 
																	 VectorImageType::Pointer field) {
	VectorImageType::Pointer image = VectorImageType::New();
	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = std::sqrt((float)pixelCount);  // size along X
	size[1]  = std::sqrt((float)pixelCount);  // size along Y
	size[2]	 = 1;
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	if(pixelCount == 1) {
		VectorImageType::IndexType index;
		index[0] = 0; index[1] = 0; index[2] = 0;
		VectorImageType::PixelType value;
		value[0] = 0; value[1] = 0; value[2] = 0;
		image->SetPixel(index, value);
	}
	else {
		itk::ImageRegionIterator<VectorImageType> iter(image, image->GetRequestedRegion());
		for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
			VectorImageType::IndexType index = iter.GetIndex();
			float smaller = field->GetRequestedRegion().GetSize()[0];
			float bigger = image->GetRequestedRegion().GetSize()[0];
			float ratio = smaller / bigger;
			int factor[3];
			factor[0] = index[0] % 2 + 1;
			factor[1] = index[1] % 2 + 1;
			factor[2] = index[2] % 2 + 1;
			index[0] = index[0] * ratio;
			index[1] = index[1] * ratio;
			index[2] = index[2] * ratio;
			VectorImageType::PixelType value = (field->GetPixel(index) * 2) / 3;
			value[0] = value[0] * factor[0];
			value[1] = value[1] * factor[1];
			value[2] = value[2] * factor[2];
			iter.Set(value);
		}
	}

	return image;
}

NLabelImageType::Pointer Utils::combineLabels(int pixelCount, GreyImageType::Pointer label) {
	NLabelImageType::Pointer image = NLabelImageType::New();
	// The image region should be initialized
	NLabelImageType::IndexType start;
	NLabelImageType::SizeType  size;
	size[0]  = std::sqrt((float)pixelCount);  // size along X
	size[1]  = std::sqrt((float)pixelCount);  // size along Y
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y

	NLabelImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	itk::ImageRegionIterator<NLabelImageType> iterImage(image, image->GetRequestedRegion());
	for(iterImage.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage) {
		NLabelImageType::PixelType value;
		for(unsigned int i = 0; i < TissueTypesMax; ++i)
			value[i] = 0;
		iterImage.Set(value);
	}

	itk::ImageRegionIterator<GreyImageType> iterLabel(label, label->GetRequestedRegion());
	for(iterLabel.GoToBegin(); !iterLabel.IsAtEnd(); ++iterLabel) {
		GreyImageType::IndexType index = iterLabel.GetIndex();
		NLabelImageType::IndexType indexV;
		float bigger = label->GetRequestedRegion().GetSize()[0];
		float smaller = image->GetRequestedRegion().GetSize()[0];
		float ratio = smaller / bigger;
		indexV[0] = index[0] * ratio;
		indexV[1] = index[1] * ratio;
		NLabelImageType::PixelType value = image->GetPixel(indexV);
		++value[(int)iterLabel.Get()];
		image->SetPixel(indexV, value);
	}

	for(iterImage.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage) {
		NLabelImageType::PixelType value = iterImage.Get();
		float total = 0;
		for(unsigned int i = 0; i < TissueTypesMax; ++i)
			total += value[i];
		float ratio = 1.0/total;
		for(unsigned int i = 0; i < TissueTypesMax; ++i)
			value[i] = value[i] * ratio;
		iterImage.Set(value);
	}

	/*
	// debug
	NeighborhoodIterator<VectorImageType>::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;
	radius[2] = 0;
	NeighborhoodIterator<VectorImageType,ZeroFluxNeumannBoundaryCondition<VectorImageType> > 
		labelIter(radius, image, image->GetRequestedRegion());
	for(labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter)
		std::cout << "[label] index " << labelIter.GetIndex() << " values " << labelIter.GetPixel(3) <<
			" " << labelIter.GetPixel(4) << " " << labelIter.GetPixel(5) << std::endl;
	*/

	return image;
}

FloatImageType::Pointer Utils::extractYM(NLabelImageType::Pointer image, vector<MechanicalProperties*>& mech) {
	FloatImageType::Pointer output = FloatImageType::New();
	output->SetRegions( image->GetBufferedRegion() );
	output->Allocate();
	itk::ImageRegionIterator<NLabelImageType> iterImage(image, image->GetRequestedRegion());
	itk::ImageRegionIterator<FloatImageType> iterOutput(output, output->GetRequestedRegion());
	for(iterImage.GoToBegin(), iterOutput.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage, ++iterOutput) {
		FloatImageType::PixelType value = 0;
		for(unsigned int i = 0; i < mech.size(); ++i)
			value += iterImage.Get()[i]*mech[i]->getYoungsModulus();
		iterOutput.Set(value);
	}
	return output;
}

VectorImageType::Pointer Utils::combineLabelsY(int pixelCount, VectorImageType::Pointer label) {
	VectorImageType::Pointer image = VectorImageType::New();
	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = std::sqrt((float)pixelCount);  // size along X
	size[1]  = std::sqrt((float)pixelCount);  // size along Y
	size[2]	 = 1;
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	itk::ImageRegionIterator<VectorImageType> iterImage(image, image->GetRequestedRegion());
	for(iterImage.GoToBegin(); !iterImage.IsAtEnd(); ++iterImage) {
		VectorImageType::PixelType value;
		value[0] = 0;
		value[1] = 0;
		value[2] = 0;
		iterImage.Set(value);
	}

	itk::ImageRegionIterator<VectorImageType> iterLabel(label, label->GetRequestedRegion());
	for(iterLabel.GoToBegin(); !iterLabel.IsAtEnd(); ++iterLabel) {
		ImageType::IndexType index = iterLabel.GetIndex();
		float bigger = label->GetRequestedRegion().GetSize()[0];
		float smaller = image->GetRequestedRegion().GetSize()[0];
		float ratio = smaller / bigger;
		index[0] = index[0] * ratio;
		index[1] = index[1] * ratio;
		index[2] = index[2] * ratio;
		VectorImageType::PixelType value = image->GetPixel(index);
		image->SetPixel(index, value);
	}

	return image;
}

VectorImageType::Pointer Utils::createBoundaryCondtionInitialDVF(int imageSize[], 
																 const float maxDeform, 
																 const int direction) {
	VectorImageType::Pointer image = VectorImageType::New();
	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();
	
	itk::ImageRegionIterator<VectorImageType> iter(image, image->GetRequestedRegion());
	switch(direction) {
		case 0:
			for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
				VectorImageType::PixelType value;
				value[0] = 0; value[1] = 0; value[2] = 0;
				VectorImageType::IndexType index = iter.GetIndex();
				if(index[0]==300) value[0] = -maxDeform;
				iter.Set(value);
			}
			break;
	}

	return image;	
}

ImageType::Pointer Utils::oneDirectionDVF(const int imageSize[], const float maxDeform) {
	ImageType::Pointer image = ImageType::New();
	// The image region should be initialized
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	// f(x,y) = (x - x0)^2 + (y - y0)^2 + z0
	//ImageType::PixelType normalise = size[0]*size[0];
	ImageType::PixelType normalise = size[0]-1;
	for(int k = 0; k < size[2]; k++) {
		for(int j = 0; j < size[1]; j++) {
			for(int i = 0; i < size[0]; i++) {
				ImageType::IndexType index;
				index[0] = i; index[1] = j; index[2] = k;
				//ImageType::PixelType elastic = i*i;
				ImageType::PixelType elastic = i;
				elastic = elastic / normalise;
				elastic = elastic * maxDeform;
				elastic *= -1;
				image->SetPixel(index, elastic);
			}
		}
	}

	return image;
}

ImageType::Pointer Utils::twoDirectionDVF(const int imageSize[], const float maxDeform) {
	ImageType::Pointer image = ImageType::New();
	// The image region should be initialized
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	// f(x,y) = (x - x0)^2 + (y - y0)^2 + z0
	//ImageType::PixelType normalise = size[1]*size[1];
	ImageType::PixelType normalise = (size[1]-1)/2;
	for(int k = 0; k < size[2]; k++) {
		for(int j = 0; j < size[1]; j++) {
			for(int i = 0; i < size[0]; i++) {
				ImageType::IndexType index;
				index[0] = i; index[1] = j; index[2] = k;
				//ImageType::PixelType elastic = (j-size[1]/2)*(j-size[1]/2);
				ImageType::PixelType half = (size[1]-1)/2;
				ImageType::PixelType elastic = j-half;
				elastic = elastic / normalise;
				elastic = elastic * maxDeform;
				image->SetPixel(index, elastic);
			}
		}
	}

	return image;
}

ImageType::Pointer Utils::createLabelsOneTissue(int imageSize[]) {
	ImageType::Pointer image = ImageType::New();
	// The image region should be initialized
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	for(int k = 0; k < size[2]; k++) {
		for(int j = 0; j < size[1]; j++) {
			for(int i = 0; i < size[0]; i++) {
				ImageType::IndexType index;
				index[0] = i; index[1] = j; index[2] = k;
				ImageType::PixelType currentPixel;

				if(j >= 50 && i >= 0 && j <= 250 && i <= 200) currentPixel = 1;
				// air
				else currentPixel = 0;

				image->SetPixel(index, currentPixel);
			}
		}
	}

	return image;
}

ImageType::Pointer Utils::createLabelsTwoTissue(int imageSize[]) {
	ImageType::Pointer image = ImageType::New();
	// The image region should be initialized
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	for(int k = 0; k < size[2]; k++) {
		for(int j = 0; j < size[1]; j++) {
			for(int i = 0; i < size[0]; i++) {
				ImageType::IndexType index;
				index[0] = i; index[1] = j; index[2] = k;
				ImageType::PixelType currentPixel;

				if(j >= 50 && i >= 0 && j <= 150 && i <= 100) currentPixel = 2;
				else if(j >= 50 && i >= 0 && j <= 250 && i <= 200) currentPixel = 1;
				// air
				else currentPixel = 0;

				image->SetPixel(index, currentPixel);
			}
		}
	}

	return image;
}

ImageType::Pointer Utils::createLabelsTwoTissueHorizontal(int imageSize[]) {
	ImageType::Pointer image = ImageType::New();
	// The image region should be initialized
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	for(int k = 0; k < size[2]; k++) {
		for(int j = 0; j < size[1]; j++) {
			for(int i = 0; i < size[0]; i++) {
				ImageType::IndexType index;
				index[0] = i; index[1] = j; index[2] = k;
				ImageType::PixelType currentPixel;

				if(j >= 50 && i >= 0 && j <= 150 && i <= 200) currentPixel = 1;
				else if(j >= 50 && i >= 0 && j <= 250 && i <= 200) currentPixel = 2;
				// air
				else currentPixel = 0;

				image->SetPixel(index, currentPixel);
			}
		}
	}

	return image;
}

ImageType::Pointer Utils::createLabelsTwoTissueSimple(int imageSize[]) {
	ImageType::Pointer image = ImageType::New();
	// The image region should be initialized
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	for(int k = 0; k < size[2]; k++) {
		for(int j = 0; j < size[1]; j++) {
			for(int i = 0; i < size[0]; i++) {
				ImageType::IndexType index;
				index[0] = i; index[1] = j; index[2] = k;
				ImageType::PixelType currentPixel;

				if(j >= 50 && i >= 0 && j <= 250 && i <= 100) currentPixel = 2;
				else if(j >= 50 && i >= 0 && j <= 250 && i <= 200) currentPixel = 1;
				// air
				else currentPixel = 0;

				image->SetPixel(index, currentPixel);
			}
		}
	}

	return image;
}

ImageType::Pointer Utils::constantDVF(int imageSize[], ImageType::PixelType constant) {
	ImageType::Pointer image = ImageType::New();
	// The image region should be initialized
	ImageType::IndexType start;
	ImageType::SizeType  size;
	size[0]  = imageSize[0];  // size along X
	size[1]  = imageSize[1];  // size along Y
	size[2]	 = imageSize[2];
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y
	start[2] =   0;

	ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();
	image->FillBuffer(constant);
	return image;
}

vtkImageData* Utils::magnitude(vtkImageData* image) {
	int dim[3];
	image->GetDimensions(dim);

	vtkImageData *magnitude=vtkImageData::New();
	magnitude->SetDimensions(dim);
	magnitude->SetOrigin( 0.0, 0.0, 0.0);
	magnitude->SetSpacing( 1.0, 1.0, 1.0);
	magnitude->SetScalarTypeToFloat();
	magnitude->SetNumberOfScalarComponents(1);
	magnitude->AllocateScalars();

	for(int i = 0; i < dim[0]; i++) {
		for(int j = 0; j < dim[1]; j++) {
			for(int k = 0; k < dim[2]; k++) {
				// get scalar pointer to current pixel
				float *x = (float*) image->GetScalarPointer( i , j , k );
				float *y = x++;
				float *z = x++;
				float *current = (float*) magnitude->GetScalarPointer( i , j , k );
				*current = sqrt((*x * *x) + (*y * *y) + (*z * *z));
			}
		}
	}

	return magnitude;
}

VectorImageType::Pointer Utils::sampleDVF(VectorImageType::Pointer image) {
	VectorImageType::Pointer sampledImage = VectorImageType::New();
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = 21; size[1]  = 21; size[2] = 1;
	start[0] = 0; start[1] = 0; start[2] = 0;
	VectorImageType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);
	sampledImage->SetRegions( region );
	sampledImage->Allocate();

	for(unsigned int j = 50, sJ = 0; j <= 250; j+=10, ++sJ) {
		for(unsigned int i = 0, sI = 0; i <= 200; i+=10, ++sI) {
			VectorImageType::IndexType pixelIndex;
			pixelIndex[0] = i; pixelIndex[1] = j; pixelIndex[2] = 0;
			VectorImageType::IndexType pixelIndexS;
			pixelIndexS[0] = sI; pixelIndexS[1] = sJ; pixelIndexS[2] = 0;
			VectorImageType::PixelType pixelValue = image->GetPixel( pixelIndex );
			sampledImage->SetPixel( pixelIndexS, pixelValue );
		}
	}

	return sampledImage;
}

void Utils::rootMeanSquareDeviation(const Vector3ImageType::Pointer image0, 
									  const Vector3ImageType::Pointer image1, 
									  const BinaryImageType::Pointer binaryImage,
									  float &rmsd) {
	rmsd = 0;
	float diff = 0;
	Vector3ImageType::PixelType v;
	itk::ImageRegionIterator<Vector3ImageType> iterImage0(image0, image0->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> iterImage1(image1, image1->GetRequestedRegion());
	itk::ImageRegionIterator<BinaryImageType> iterBinary(binaryImage, binaryImage->GetRequestedRegion());
	unsigned int counter = 0;
	for(iterImage0.GoToBegin(), iterImage1.GoToBegin(), iterBinary.GoToBegin(); 
			!iterImage0.IsAtEnd(); ++iterImage0, ++iterImage1, ++iterBinary, ++counter) {
		if(iterBinary.Get() == 1) {
			v = iterImage0.Get() - iterImage1.Get();
			diff = std::sqrt(v[0]*v[0] + v[1]*v[1]);
			rmsd += diff*diff;
		}
	}
	rmsd = std::sqrt(rmsd / counter);
}

void Utils::printDifference(const Vector3ImageType::Pointer image0, 
							const Vector3ImageType::Pointer image1) {
	itk::ImageRegionIterator<Vector3ImageType> iterImage0(image0, image0->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> iterImage1(image1, image1->GetRequestedRegion());
	for(; !iterImage0.IsAtEnd(); ++iterImage0, ++iterImage1) {
		for(unsigned int i = 0; i < image0->GetImageDimension(); i++) {
			double diff = iterImage0.Get()[i] - iterImage1.Get()[i];
			if(diff != 0)
				cout << iterImage0.GetIndex() << " " << diff << endl;
		}
	}
}

Vector3ImageType::Pointer Utils::readAbaqusFile(const char* fileName) {
	// Open the input file.
	ifstream fin(fileName);
	if(!fin) {
		cerr << "File not found" << endl;
		return 0;
    }

	Vector3ImageType::Pointer image = Vector3ImageType::New();

	// The image region should be initialized
	Vector3ImageType::IndexType start;
	Vector3ImageType::SizeType  size;
	size[0]  = 61;	size[1]  = 61; size[2] = 1;
	start[0] = 0; start[1] = 0; start[2] = 0;

	Vector3ImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	float translation[3] = {50,250,0};
	float shrink[3] = {0.2,0.2,1};
	float scale[3] = {-1,1,1};

	for(unsigned int i = 0; i < size[0]*size[1]*size[2]; i++) {
		Vector3ImageType::PixelType value;
		std::string comment;
		fin >> comment;
		float x[5];
		fin >> x[0] >> x[1] >> x[2] >> x[3] >> x[4];

		Vector3ImageType::IndexType index;
		for(unsigned int j = 0; j < 3; j++)
			index[j] = (x[j]+translation[j])*shrink[j];

		value[0] = x[3]*scale[0];
		value[1] = x[4]*scale[1];
		value[2] = 0*scale[2];

		//std::cout << "index: " << index << ", value: " << value << std::endl;
		image->SetPixel(index, value);
	}

	return image;
}

void Utils::splitImage(Vector3ImageType::Pointer image, FloatImageType::Pointer& image1,
											 FloatImageType::Pointer& image2) {
	image1 = FloatImageType::New();
	image2 = FloatImageType::New();

	Vector3ImageType::RegionType  region3D = image->GetBufferedRegion();
	Vector3ImageType::IndexType   index3D = region3D.GetIndex();;
	Vector3ImageType::SizeType    size3D = region3D.GetSize();
	
	FloatImageType::RegionType  region;
	FloatImageType::IndexType   index;
	FloatImageType::SizeType    size;

	index[0] = index3D[0];
	index[1] = index3D[1];
	size[0]  = size3D[0];
	size[1]  = size3D[1];
	region.SetSize( size );
	region.SetIndex( index );

	image1->SetRegions( region );
	image2->SetRegions( region );
	image1->Allocate();
	image2->Allocate();

	itk::ImageRegionIterator<FloatImageType> floatIter1(image1, image1->GetRequestedRegion());
	itk::ImageRegionIterator<FloatImageType> floatIter2(image2, image2->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> imageIter(image, image->GetRequestedRegion());
	for(imageIter.GoToBegin(), floatIter1.GoToBegin(), floatIter2.GoToBegin(); 
		!imageIter.IsAtEnd(); ++imageIter, ++floatIter1, ++floatIter2) {
			floatIter1.Set(imageIter.Get()[0]);
			floatIter2.Set(imageIter.Get()[1]);
	}
}

void Utils::splitImage(VectorImageType::Pointer image, FloatImageType::Pointer& image1,
											 FloatImageType::Pointer& image2) {
	image1 = FloatImageType::New();
	image2 = FloatImageType::New();
	image1->SetRegions(image->GetBufferedRegion());
	image2->SetRegions(image->GetBufferedRegion());
	image1->Allocate();
	image2->Allocate();

	itk::ImageRegionIterator<FloatImageType> floatIter1(image1, image1->GetRequestedRegion());
	itk::ImageRegionIterator<FloatImageType> floatIter2(image2, image2->GetRequestedRegion());
	itk::ImageRegionIterator<VectorImageType> imageIter(image, image->GetRequestedRegion());
	for(imageIter.GoToBegin(), floatIter1.GoToBegin(), floatIter2.GoToBegin(); 
		!imageIter.IsAtEnd(); ++imageIter, ++floatIter1, ++floatIter2) {
			floatIter1.Set(imageIter.Get()[0]);
			floatIter2.Set(imageIter.Get()[1]);
	}
}

void Utils::combineImage(FloatImageType::Pointer image1, FloatImageType::Pointer image2,
						 Vector3ImageType::Pointer& image) {
	image = Vector3ImageType::New();

	FloatImageType::RegionType  region = image1->GetBufferedRegion();
	FloatImageType::IndexType   index = region.GetIndex();
	FloatImageType::SizeType    size = region.GetSize();

	Vector3ImageType::RegionType  region3D;
	Vector3ImageType::IndexType   index3D;
	Vector3ImageType::SizeType    size3D;
	
	index3D[0] = index[0];
	index3D[1] = index[1];
	index3D[2] = 0;
	size3D[0]  = size[0];
	size3D[1]  = size[1];
	size3D[2]  = 1;
	region3D.SetSize( size3D );
	region3D.SetIndex( index3D );

	image->SetRegions( region3D );
	image->Allocate();

	itk::ImageRegionIterator<FloatImageType> floatIter1(image1, image1->GetRequestedRegion());
	itk::ImageRegionIterator<FloatImageType> floatIter2(image2, image2->GetRequestedRegion());
	itk::ImageRegionIterator<Vector3ImageType> imageIter(image, image->GetRequestedRegion());
	for(imageIter.GoToBegin(), floatIter1.GoToBegin(), floatIter2.GoToBegin(); 
		!imageIter.IsAtEnd(); ++imageIter, ++floatIter1, ++floatIter2) {
			Vector3ImageType::PixelType value;
			value[0] = floatIter1.Get();
			value[1] = floatIter2.Get();
			value[2] = 0;
			imageIter.Set(value);
	}
}

void Utils::combineImage(FloatImageType::Pointer image1, FloatImageType::Pointer image2,
						 VectorImageType::Pointer& image) {
	image = VectorImageType::New();
	image->SetRegions( image1->GetBufferedRegion() );
	image->Allocate();

	itk::ImageRegionIterator<FloatImageType> floatIter1(image1, image1->GetRequestedRegion());
	itk::ImageRegionIterator<FloatImageType> floatIter2(image2, image2->GetRequestedRegion());
	itk::ImageRegionIterator<VectorImageType> imageIter(image, image->GetRequestedRegion());
	for(imageIter.GoToBegin(), floatIter1.GoToBegin(), floatIter2.GoToBegin(); 
		!imageIter.IsAtEnd(); ++imageIter, ++floatIter1, ++floatIter2) {
			VectorImageType::PixelType value;
			value[0] = floatIter1.Get();
			value[1] = floatIter2.Get();
			imageIter.Set(value);
	}
}

void Utils::testBoundaryCondition() {
	VectorImageType::Pointer image = VectorImageType::New();
	// The image region should be initialized
	VectorImageType::IndexType start;
	VectorImageType::SizeType  size;
	size[0]  = 1;  // size along X
	size[1]  = 1;  // size along Y
	start[0] =   0;  // first index on X
	start[1] =   0;  // first index on Y

	VectorImageType::RegionType region;
	region.SetSize( size );
	region.SetIndex( start );

	// Pixel data is allocated
	image->SetRegions( region );
	image->Allocate();

	VectorImageType::IndexType index;
	index[0] = 0; index[1] = 0;
	VectorImageType::PixelType value;
	value[0] = 0; value[1] = 0;
	image->SetPixel(index, value);

	typedef DisplacementBoundaryCondition<VectorImageType> BC;
	BC bc;
	VectorPixelType constant;
	constant[0] = 0; constant[1] = 0;
	VectorPixelType displacementRight;
	displacementRight[0] = -15; displacementRight[1] = 0;
	VectorPixelType displacementLeft;
	displacementLeft[0] = 15; displacementLeft[1] = 0;
	VectorPixelType displacementTop;
	displacementTop[0] = 0; displacementTop[1] = 10;
	VectorPixelType displacementBottom;
	displacementBottom[0] = 0; displacementBottom[1] = -10;
	bc.SetConstant(constant);
	bc.SetMaxDisplacementRight(displacementRight);
	bc.SetMaxDisplacementLeft(displacementLeft);
	bc.SetMaxDisplacementTop(displacementTop);
	bc.SetMaxDisplacementBottom(displacementBottom);

	NeighborhoodIterator<VectorImageType>::RadiusType radius;
	radius[0] = 1;
	radius[1] = 1;

	NeighborhoodIterator<VectorImageType,BC> 
		bit(radius, image, image->GetRequestedRegion());
	bit.OverrideBoundaryCondition(&bc);

	for(bit.GoToBegin(); !bit.IsAtEnd(); ++bit) {
		std::cout << "bc top " << bit.GetPixel(1) << std::endl;
		std::cout << "bc left " << bit.GetPixel(3) << std::endl;
		std::cout << "bc right " << bit.GetPixel(5) << std::endl;
		std::cout << "bc bottom " << bit.GetPixel(7) << std::endl;
	}
}

GreyImageType::Pointer Utils::resampling(int numberOfPixels, GreyImageType::Pointer input) {
	typedef itk::ResampleImageFilter<GreyImageType,GreyImageType> ResampleFilterType;
	ResampleFilterType::Pointer filter = ResampleFilterType::New();
	typedef itk::AffineTransform< double, 2 > TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->Scale((input->GetBufferedRegion().GetSize()[0]-1.0)/(std::sqrt((float)numberOfPixels)-1.0));
	filter->SetTransform( transform );
	typedef itk::NearestNeighborInterpolateImageFunction<GreyImageType, double > InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	filter->SetInterpolator( interpolator );
	double spacing[ Dimension ];
	spacing[0] = 1.0; spacing[1] = 1.0;
	filter->SetOutputSpacing( spacing );
	double origin[ Dimension ];
	origin[0] = 0.0; origin[1] = 0.0;
	filter->SetOutputOrigin( origin );
	GreyImageType::SizeType size;
	size[0] = std::sqrt((float)numberOfPixels); size[1] = std::sqrt((float)numberOfPixels);
	filter->SetSize( size );
	filter->SetInput( input );
	filter->Update();
	return filter->GetOutput();
}

VectorImageType::Pointer Utils::resampling(int pixels, VectorImageType::Pointer input) {
	if(pixels == 1) {
		VectorImageType::Pointer image = VectorImageType::New();
		VectorImageType::IndexType start;
		VectorImageType::SizeType  size;
		size[0] = 1; size[1] = 1;
		start[0] = 0; start[1] = 0;
		VectorImageType::RegionType region;
		region.SetSize( size );
		region.SetIndex( start );
		image->SetRegions( region );
		image->Allocate();
		VectorImageType::IndexType index;
		index[0] = 0; index[1] = 0;
		VectorImageType::PixelType value;
		value[0] = 0; value[1] = 0;
		image->SetPixel(index, value);
		return image;
	}

	// resample
	FloatImageType::Pointer single[2];
	FloatImageType::Pointer singleOut[2];
	Utils::splitImage(input, single[0], single[1]);
	for(int i = 0; i < 2; ++i) {
		typedef itk::ResampleImageFilter<FloatImageType,FloatImageType> ResampleFilterType;
		ResampleFilterType::Pointer filter = ResampleFilterType::New();
		typedef itk::AffineTransform< double, 2 > TransformType;
		TransformType::Pointer transform = TransformType::New();
		transform->Scale((input->GetBufferedRegion().GetSize()[0]-1.0)/(std::sqrt((float)pixels)-1.0));
		filter->SetTransform( transform );
		typedef itk::LinearInterpolateImageFunction<FloatImageType, double > InterpolatorType;
		//typedef itk::NearestNeighborInterpolateImageFunction<FloatImageType, double > InterpolatorType;
		InterpolatorType::Pointer interpolator = InterpolatorType::New();
		filter->SetInterpolator( interpolator );
		double spacing[ Dimension ];
		spacing[0] = 1.0; spacing[1] = 1.0;
		filter->SetOutputSpacing( spacing );
		double origin[ Dimension ];
		origin[0] = 0.0; origin[1] = 0.0;
		filter->SetOutputOrigin( origin );
		FloatImageType::SizeType size;
		size[0] = std::sqrt((float)pixels); size[1] = std::sqrt((float)pixels);
		filter->SetSize( size );
		filter->SetInput( single[i] );
		filter->Update();
		singleOut[i] = filter->GetOutput();
	}
	VectorImageType::Pointer combined;
	Utils::combineImage(singleOut[0], singleOut[1], combined);
	return combined;
}
