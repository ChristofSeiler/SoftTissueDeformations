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

#include "vtkAbaqusPointsReader.h"

#include "vtkCellArray.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

vtkCxxRevisionMacro(vtkAbaqusPointsReader, "$Revision: 1.0.0.0 $");
vtkStandardNewMacro(vtkAbaqusPointsReader);

//----------------------------------------------------------------------------
vtkAbaqusPointsReader::vtkAbaqusPointsReader()
{
  this->FileName = 0;
  this->SetNumberOfInputPorts(0);
}

//----------------------------------------------------------------------------
vtkAbaqusPointsReader::~vtkAbaqusPointsReader()
{
  this->SetFileName(0);
}

//----------------------------------------------------------------------------
void vtkAbaqusPointsReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n";

}

//----------------------------------------------------------------------------
int vtkAbaqusPointsReader::RequestData(vtkInformation*,
                                       vtkInformationVector**,
                                       vtkInformationVector* outputVector)
{
  // Make sure we have a file to read.
  if(!this->FileName)
    {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }

  // Open the input file.
  ifstream fin(this->FileName);
  if(!fin)
    {
    vtkErrorMacro("Error opening file " << this->FileName);
    return 0;
    }

  // Assign points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
  scalars->SetNumberOfComponents(1);
  vtkSmartPointer<vtkFloatArray> vectors = vtkSmartPointer<vtkFloatArray>::New();
  vectors->SetNumberOfComponents(3);

  int res[3] = { 21, 21, 1 };
  int noOfTubles = res[0]*res[1]*res[2];
  scalars->SetNumberOfTuples(noOfTubles);
  vectors->SetNumberOfTuples(noOfTubles);
  points->Allocate(noOfTubles);

  int id = 0;

  for ( int k=0; k<res[2]; ++k) {
	  for ( int j=0; j<res[1]; ++j) {
		  for ( int i=0; i<res[0]; ++i) {
			  float x[3];
			  fin >> x[0] >> x[1] >> x[2];
			  float v[3];
			  v[0] = x[1]; v[1] = x[2]; v[2] = 0;
			  // left bottom corner is origin in vtk
			  points->InsertPoint(id, i, res[1]-1-j, k);
			  float magnitude[1] = { sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) };
			  scalars->InsertTuple(id, magnitude);
			  vectors->InsertTuple(id, v);
			  id++;
		  }
	  }
  }

  // Store the points and cells in the output data object.
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  // Assign scalars
  output->GetPointData()->SetScalars(scalars);
  // Assign vectors
  output->GetPointData()->SetVectors(vectors);

  return 1;
}
