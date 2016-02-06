PROJECT(BCPreperation)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

INCLUDE_DIRECTORIES(${BCPreperation_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES($ENV{INSIGHT_APPLICATIONS}/Auxiliary/vtk)
INCLUDE_DIRECTORIES($ENV{CG_INC_PATH})
INCLUDE_DIRECTORIES(${BCPreperation_SOURCE_DIR}/../DVFRegularization/includes)
INCLUDE_DIRECTORIES(${BCPreperation_SOURCE_DIR}/../DVFRegularization/src)

ADD_EXECUTABLE(
BCPreperation  
${BCPreperation_SOURCE_DIR}/src/Main.cxx

)

TARGET_LINK_LIBRARIES(BCPreperation ITKBasicFilters ITKCommon ITKIO vtkRendering 
vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon vtkWidgets itkStatistics
$ENV{CG_LIB_PATH}/cg.lib $ENV{CG_LIB_PATH}/cgGL.lib glew32.lib)
