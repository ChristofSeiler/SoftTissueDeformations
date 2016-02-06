PROJECT(ParameterEstimation)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

INCLUDE_DIRECTORIES(${ParameterEstimation_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES(E:/libraries/itk/InsightApplications-3.2.0/Auxiliary/vtk)
INCLUDE_DIRECTORIES(${ParameterEstimation_SOURCE_DIR}/../DVFRegularization/includes)
INCLUDE_DIRECTORIES(${ParameterEstimation_SOURCE_DIR}/../DVFRegularization/src)

ADD_EXECUTABLE(
ParameterEstimation  
${ParameterEstimation_SOURCE_DIR}/includes/MetropolisSamplerMinimizer.h
${ParameterEstimation_SOURCE_DIR}/includes/RegularizationCostFunction.h
${ParameterEstimation_SOURCE_DIR}/src/Main.cxx
${ParameterEstimation_SOURCE_DIR}/src/MetropolisSamplerMinimizer.cxx
${ParameterEstimation_SOURCE_DIR}/src/RegularizationCostFunction.cxx

${ParameterEstimation_SOURCE_DIR}/../DVFRegularization/src/Utils.cxx
${ParameterEstimation_SOURCE_DIR}/../DVFRegularization/src/MechanicalProperties.cxx
${ParameterEstimation_SOURCE_DIR}/../DVFRegularization/src/vtkAbaqusPointsReader.cxx

)

TARGET_LINK_LIBRARIES(ParameterEstimation ITKBasicFilters ITKCommon ITKIO vtkRendering vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon vtkWidgets itkStatistics)
