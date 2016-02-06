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

#ifndef __vtkAbaqusPointsReader_h
#define __vtkAbaqusPointsReader_h

#include "vtkPolyDataAlgorithm.h"

class VTK_IO_EXPORT vtkAbaqusPointsReader : public vtkPolyDataAlgorithm
{
public:
  static vtkAbaqusPointsReader* New();
  vtkTypeRevisionMacro(vtkAbaqusPointsReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get the name of the file from which to read points.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

protected:
  vtkAbaqusPointsReader();
  ~vtkAbaqusPointsReader();

  char* FileName;

  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
private:
  vtkAbaqusPointsReader(const vtkAbaqusPointsReader&);  // Not implemented.
  void operator=(const vtkAbaqusPointsReader&);  // Not implemented.
};

#endif
