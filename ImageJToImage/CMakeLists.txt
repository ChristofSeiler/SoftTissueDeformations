PROJECT(ImageJToImage)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${USE_ITK_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(VTK)
IF(VTK_FOUND)
INCLUDE(${USE_VTK_FILE})
ENDIF(VTK_FOUND)

INCLUDE_DIRECTORIES(${ImageJToImage_SOURCE_DIR}/includes)
INCLUDE_DIRECTORIES($ENV{INSIGHT_APPLICATIONS}/Auxiliary/vtk)
INCLUDE_DIRECTORIES($ENV{CG_INC_PATH})
INCLUDE_DIRECTORIES(${ImageJToImage_SOURCE_DIR}/../DVFRegularization/includes)
INCLUDE_DIRECTORIES(${ImageJToImage_SOURCE_DIR}/../DVFRegularization/src)

ADD_EXECUTABLE(
ImageJToImage  
${ImageJToImage_SOURCE_DIR}/src/Main.cxx

)

TARGET_LINK_LIBRARIES(ImageJToImage ITKBasicFilters ITKCommon ITKIO vtkRendering 
vtkGraphics vtkHybrid vtkImaging vtkIO vtkFiltering vtkCommon vtkWidgets itkStatistics
$ENV{CG_LIB_PATH}/cg.lib $ENV{CG_LIB_PATH}/cgGL.lib $ENV{CG_LIB_PATH}/glew32.lib)
