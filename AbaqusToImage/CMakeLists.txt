PROJECT(AbaqusToImage)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

INCLUDE_DIRECTORIES(${AbaqusToImage_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES($ENV{INSIGHT_APPLICATIONS}/Auxiliary/vtk)
INCLUDE_DIRECTORIES($ENV{CG_INC_PATH})
INCLUDE_DIRECTORIES(${AbaqusToImage_SOURCE_DIR}/../DVFRegularization/includes)
INCLUDE_DIRECTORIES(${AbaqusToImage_SOURCE_DIR}/../DVFRegularization/src)

ADD_EXECUTABLE(
AbaqusToImage  
${AbaqusToImage_SOURCE_DIR}/src/Main.cxx
${AbaqusToImage_SOURCE_DIR}/../DVFRegularization/src/Utils.cxx

)

TARGET_LINK_LIBRARIES(AbaqusToImage ITKBasicFilters ITKCommon ITKIO vtkRendering 
vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon vtkWidgets itkStatistics
$ENV{CG_LIB_PATH}/cg.lib $ENV{CG_LIB_PATH}/cgGL.lib glew32.lib)
