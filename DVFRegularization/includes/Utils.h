/*
 * Displacement Vector Field Regularization for Modelling of Soft Tissue Deformations
 * Master Thesis
 * Master of Science in Biomedical Engineering
 * University of Bern
 * \author	Christof Seiler
 * \date	September 07 - May 08
 * \brief	Utils for general use.
 */

#pragma once

#include "Definitions.h"
#include <vector>

// foward declaration
class vtkImageData;
class vtkPolyData;
class vtkPolyDataMapper;
namespace dvfRegularization {
	class MechanicalProperties;
}

using namespace std;

namespace dvfRegularization {

	struct Count {
		Count() 
			: horizontalLine(NULL), verticalLine(NULL) {}
		Count(std::vector<int>* hl, std::vector<int>* vl) 
			: horizontalLine(hl), verticalLine(vl) {}
		~Count() {
			if(horizontalLine != NULL)
				delete horizontalLine;
			if(verticalLine != NULL)
				delete verticalLine;
		}
		std::vector<int>* horizontalLine;
		std::vector<int>* verticalLine;
	};

	class Utils
	{
	public:
		static void visualizeImage(vtkImageData *im);
		static void visualizeMultipleImages(vector<vtkImageData*> &images);
		static void visualizeMultipleImages3D(vector<vtkImageData*> &images);
		static void visualizeImage3D(vtkImageData *im);
		static void visualizeDvfArrow(vtkPolyData *image, const float factor);
		static void compoundImage(vtkImageData *dvfx,vtkImageData *dvfy,vtkImageData *dvfz,vtkImageData *DVF);
		static void deformImage(vtkImageData *inputImage,vtkImageData *DVF,vtkImageData *outputImage);
		static vtkImageData* scalarToColor(vtkImageData* image);
		static vtkPolyData* convertImageToPoly(vtkImageData* image);
		static vtkPolyData* convertImageToPoly(const VectorImageType::Pointer image);
		static VectorImageType::Pointer convertPolyToITK(vtkPolyData* poly);
		static VectorImageType::Pointer createVectorImage(ImageType::Pointer image[3]);
		static vtkImageData* convertSampledVectorImage(VectorImageType::Pointer image);
		static vtkImageData* convertVectorImage(VectorImageType::Pointer image);
		static void createGridLineImage(vtkImageData* image, int gridSize[]);
		static ImageType::Pointer createInitalElasticDVF(int imageSize[], float maxDeform);
		static VectorImageType::Pointer createComplicatedInitialGuessDVF(const int imageSize[], 
			const float maxDeform, 
			const ImageType::Pointer labels, 
			const std::vector<MechanicalProperties*> &mechanicalProperties);
		static VectorImageType::Pointer extrapolateNearestNeighborInitialDVF(int pixelCount, 
			VectorImageType::Pointer field);
		static VectorImageType::Pointer extrapolateLinearInitialDVF(int pixelCount, 
			VectorImageType::Pointer field);
		static NLabelImageType::Pointer combineLabels(int pixelCount, GreyImageType::Pointer label);
		static FloatImageType::Pointer extractYM(NLabelImageType::Pointer image, 
			vector<MechanicalProperties*>& mech);
		static VectorImageType::Pointer combineLabelsY(int pixelCount, VectorImageType::Pointer label);
		static VectorImageType::Pointer createEmtpyDVF(int imageSize[]);
		static VectorImageType::Pointer createBoundaryCondtionInitialDVF(int imageSize[], 
			const float maxDeform, const int direction);
		static ImageType::Pointer oneDirectionDVF(const int imageSize[], const float maxDeform);
		static ImageType::Pointer twoDirectionDVF(const int imageSize[], const float maxDeform);
		static ImageType::Pointer createLabelsOneTissue(int imageSize[]);
		static ImageType::Pointer createLabelsTwoTissue(int imageSize[]);
		static ImageType::Pointer createLabelsTwoTissueHorizontal(int imageSize[]);
		static ImageType::Pointer createLabelsTwoTissueSimple(int imageSize[]);
		static ImageType::Pointer constantDVF(int imageSize[], ImageType::PixelType constant);
		static vtkImageData* magnitude(vtkImageData* image);
		static VectorImageType::Pointer sampleDVF(VectorImageType::Pointer image);
		static void rootMeanSquareDeviation(const Vector3ImageType::Pointer image0, 
			const Vector3ImageType::Pointer image1, const BinaryImageType::Pointer binaryImage,
			float &rmsd);
		static void printDifference(const Vector3ImageType::Pointer image0, 
			const Vector3ImageType::Pointer image1);
		static Vector3ImageType::Pointer readAbaqusFile(const char* fileName);
		static void testBoundaryCondition();
		static void splitImage(Vector3ImageType::Pointer image, FloatImageType::Pointer& image1,
			FloatImageType::Pointer& image2);
		static void splitImage(VectorImageType::Pointer image, FloatImageType::Pointer& image1,
			FloatImageType::Pointer& image2);
		static void combineImage(FloatImageType::Pointer image1, FloatImageType::Pointer image2,
			Vector3ImageType::Pointer& image);
		static void combineImage(FloatImageType::Pointer image1, FloatImageType::Pointer image2,
			VectorImageType::Pointer& image);
		static VectorImageType::Pointer resampling(int pixels, VectorImageType::Pointer input);
		static GreyImageType::Pointer resampling(int numberOfPixels, GreyImageType::Pointer input);
	};

}