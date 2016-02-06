PROJECT(ImagePreperation)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

INCLUDE_DIRECTORIES(${ImagePreperation_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES($ENV{INSIGHT_APPLICATIONS}/Auxiliary/vtk)
INCLUDE_DIRECTORIES(${ImagePreperation_SOURCE_DIR}/../DVFRegularization/includes)
INCLUDE_DIRECTORIES(${ImagePreperation_SOURCE_DIR}/../DVFRegularization/src)

ADD_EXECUTABLE(
ImagePreperation  
${ImagePreperation_SOURCE_DIR}/src/Main.cxx

)

TARGET_LINK_LIBRARIES(ImagePreperation ITKBasicFilters ITKCommon ITKIO vtkRendering vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon vtkWidgets itkStatistics)
