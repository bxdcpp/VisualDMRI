project(vtkDMRI)

#-----------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------
if(NOT DEFINED BUILD_SHARED_LIBS)
  option(BUILD_SHARED_LIBS "Build with shared libraries." ON)
endif()

# --------------------------------------------------------------------------
# Configure headers
# --------------------------------------------------------------------------
set(configure_header_file vtkDMRIConfigure.h)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${configure_header_file}.in
  ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
  )

# --------------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------------

#
# Teem
#
find_package(Teem REQUIRED NO_MODULE)
include(${Teem_USE_FILE})

#
#VTK
#
if(${VTK_VERSION} VERSION_LESS "8.90.0")
  include(${VTK_USE_FILE})
endif()



# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------
set(vtkDMRI_SRCS
  vtkTeemEstimateDiffusionTensor.cxx
  vtkPolyDataTensorToColor.cxx
  vtkPolyDataColorLinesByOrientation.cxx
  vtkBSplineInterpolateImageFunction.cxx
  vtkHyperPointandArray.cxx
  vtkHyperStreamlineDTMRI.cxx
  vtkHyperStreamlineTeem.cxx
  vtkPreciseHyperStreamline.cxx
  vtkPreciseHyperStreamlinePoints.cxx
  vtkSeedTracts.cxx
  vtkTensorImplicitFunctionToFunctionSet.cxx
  vtkTractographyPointAndArray.cxx
  vtkTensorMask.cxx
  vtkTensorRotate.cxx
  vtkTRKWriter.cxx
  vtkTRKReader.cxx
  defs.h
  trkfileio.h
  )

# Helper classes
set_source_files_properties(
  vtkHyperPointandArray.cxx
  vtkTractographyPointAndArray.cxx
  WRAP_EXCLUDE
  )

# --------------------------------------------------------------------------
# Include dirs
# --------------------------------------------------------------------------
# set(include_dirs
#   ${CMAKE_CURRENT_SOURCE_DIR}
#   ${CMAKE_CURRENT_BINARY_DIR}
#   )

#include_directories(BEFORE ${include_dirs} ${vtkTeem_INCLUDE_DIRS})


# --------------------------------------------------------------------------
# Build library
# --------------------------------------------------------------------------
set(lib_name ${PROJECT_NAME})

set(srcs ${vtkDMRI_SRCS})

add_library(${lib_name} ${srcs})



target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${vtkTeem_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>

    PRIVATE
    $<BUILD_INTERFACE:${VTK_INCLUDE_DIRS}>
 
    )

message(Status ${CMAKE_CURRENT_BINARY_DIR})
set(libs
  vtkTeem
  #${Teem_LIBRARIES}
  ${VTK_LIBRARIES}
  )
target_link_libraries(${lib_name} ${libs})




#set_target_properties(${lib_name} PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_BIN_DIR}"
#    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}"
#    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}"
#    )

# Apply user-defined properties to the library target.
if(Slicer_LIBRARY_PROPERTIES)
  set_target_properties(${lib_name} PROPERTIES ${Slicer_LIBRARY_PROPERTIES})
endif()

# --------------------------------------------------------------------------
# Folder
# --------------------------------------------------------------------------
if(NOT DEFINED ${lib_name}_FOLDER)
  set(${lib_name}_FOLDER "Libs")
endif()
if(NOT "${${lib_name}_FOLDER}" STREQUAL "")
  set_target_properties(${lib_name} PROPERTIES FOLDER ${${lib_name}_FOLDER})
endif()






# --------------------------------------------------------------------------
# Export target
# --------------------------------------------------------------------------
if(NOT DEFINED ${lib_name}_EXPORT_FILE)
  set(${lib_name}_EXPORT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}Targets.cmake)
endif()
export(TARGETS ${lib_name} APPEND FILE ${${lib_name}_EXPORT_FILE})


# export(
#         TARGETS ${PROJECT_NAME}
#         FILE ${CMAKE_BINARY_DIR}/${PROJECT_NAME}Targets.cmake
#     )


# --------------------------------------------------------------------------
# Python Wrapping
# --------------------------------------------------------------------------
# if(VTK_WRAP_PYTHON)
#   include(vtkMacroKitPythonWrap)

#   SlicerMacroPythonWrapModuleVTKLibrary(
#     NAME ${lib_name}
#     SRCS ${vtkDMRI_SRCS}
#     RELATIVE_PYTHON_DIR "."
#     )

#   # Set python module logic output
#   set_target_properties(${lib_name}Python ${lib_name}PythonD PROPERTIES
#     RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_BIN_DIR}"
#     LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}"
#     ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${Slicer_QTLOADABLEMODULES_LIB_DIR}"
#     )
#   # Export target
#   export(TARGETS ${lib_name}Python ${lib_name}PythonD APPEND FILE ${${lib_name}_EXPORT_FILE})
# endif()

# --------------------------------------------------------------------------
# Set INCLUDE_DIRS variable
# --------------------------------------------------------------------------
set(${lib_name}_INCLUDE_DIRS ${Slicer_ACTUAL_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL "${lib_name} include dirs" FORCE)

# --------------------------------------------------------------------------
# Install library
# --------------------------------------------------------------------------
install(TARGETS ${lib_name}
  RUNTIME DESTINATION ${Slicer_INSTALL_QTLOADABLEMODULES_BIN_DIR} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${Slicer_INSTALL_QTLOADABLEMODULES_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${Slicer_INSTALL_QTLOADABLEMODULES_LIB_DIR} COMPONENT Development
  )

