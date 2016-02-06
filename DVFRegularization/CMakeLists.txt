PROJECT(DVFRegularization)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

INCLUDE_DIRECTORIES(${DVFRegularization_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES($ENV{INSIGHT_APPLICATIONS}/Auxiliary/vtk)
INCLUDE_DIRECTORIES($ENV{CG_INC_PATH})

ADD_EXECUTABLE(
DVFRegularization  
${DVFRegularization_SOURCE_DIR}/includes/Definitions.h
${DVFRegularization_SOURCE_DIR}/includes/Utils.h
${DVFRegularization_SOURCE_DIR}/includes/itkMRFRegularizationFilter.h
${DVFRegularization_SOURCE_DIR}/includes/itkMRFRegularizationFilterGPU.h
${DVFRegularization_SOURCE_DIR}/includes/itkAddGaussianNoiseFilter.h
${DVFRegularization_SOURCE_DIR}/includes/itkEnergyFunction.h
${DVFRegularization_SOURCE_DIR}/includes/itkTestEnergyFunctionImpl.h
${DVFRegularization_SOURCE_DIR}/includes/itkYoungPoissonEnergyFunction.h
${DVFRegularization_SOURCE_DIR}/includes/MechanicalProperties.h
${DVFRegularization_SOURCE_DIR}/includes/vtkAbaqusPointsReader.h
${DVFRegularization_SOURCE_DIR}/includes/APrioriImageInformation.h
${DVFRegularization_SOURCE_DIR}/includes/DisplacementBoundaryCondition.h
${DVFRegularization_SOURCE_DIR}/includes/itkInterpolateNextStepFilter.h
${DVFRegularization_SOURCE_DIR}/src/Main.cxx
${DVFRegularization_SOURCE_DIR}/src/Utils.cxx
${DVFRegularization_SOURCE_DIR}/src/itkMRFRegularizationFilter.txx
${DVFRegularization_SOURCE_DIR}/src/itkMRFRegularizationFilterGPU.txx
${DVFRegularization_SOURCE_DIR}/src/itkAddGaussianNoiseFilter.txx
${DVFRegularization_SOURCE_DIR}/src/itkEnergyFunction.txx
${DVFRegularization_SOURCE_DIR}/src/itkTestEnergyFunctionImpl.txx
${DVFRegularization_SOURCE_DIR}/src/itkYoungPoissonEnergyFunction.txx
${DVFRegularization_SOURCE_DIR}/src/itkInterpolateNextStepFilter.txx
${DVFRegularization_SOURCE_DIR}/src/MechanicalProperties.cxx
${DVFRegularization_SOURCE_DIR}/src/vtkAbaqusPointsReader.cxx
${DVFRegularization_SOURCE_DIR}/src/APrioriImageInformation.cxx
${DVFRegularization_SOURCE_DIR}/src/shader.cg
)

TARGET_LINK_LIBRARIES(DVFRegularization ITKBasicFilters ITKCommon ITKIO vtkRendering 
vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon vtkWidgets itkStatistics
$ENV{CG_LIB_PATH}/cg.lib $ENV{CG_LIB_PATH}/cgGL.lib glew32.lib)
