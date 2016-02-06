PROJECT(EnergyFunctionPlot)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

INCLUDE_DIRECTORIES(${EnergyFunctionPlot_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES(E:/libraries/itk/InsightApplications-3.2.0/Auxiliary/vtk)
INCLUDE_DIRECTORIES(${EnergyFunctionPlot_SOURCE_DIR}/../DVFRegularization/includes)
INCLUDE_DIRECTORIES(${EnergyFunctionPlot_SOURCE_DIR}/../DVFRegularization/src)

ADD_EXECUTABLE(
EnergyFunctionPlot  
${EnergyFunctionPlot_SOURCE_DIR}/src/Main.cxx

${EnergyFunctionPlot_SOURCE_DIR}/../DVFRegularization/src/Utils.cxx
${EnergyFunctionPlot_SOURCE_DIR}/../DVFRegularization/src/MechanicalProperties.cxx

)

TARGET_LINK_LIBRARIES(EnergyFunctionPlot ITKBasicFilters ITKCommon ITKIO vtkRendering vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon vtkWidgets itkStatistics)
