cmake_minimum_required (VERSION 3.14.0)
project(VisualDMRI VERSION 6.0.0)

set(VTK_SRC ${CMAKE_SOURCE_DIR}/../slicer-build/VTK)
set(VTK_DIR ${CMAKE_SOURCE_DIR}/../slicer-build/VTK-build)
set(ITK_SRC ${CMAKE_SOURCE_DIR}/../slicer-build/ITK)
set(ITK_DIR ${CMAKE_SOURCE_DIR}/../slicer-build/ITK-build)
set(Teem_SRC T:/Slicer-build/teem)
set(Teem_DIR T:/Slicer-build/teem-build)

option(BUILD_VisualNeurosurgery "Build VisualNeurosurgery" OFF)


set(CMAKE_BUILD_TYPE Debug CACHE STRING "Set the type of build.")
set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} CACHE STRING "" FORCE)
option(BUILD_SHARED_LIBS "Build the shared library" ON)
option(VTK_WRAP_PYTHON "VTK WRAP PYTHON" ON)
option(BUILD_TESTING "Test the project" OFF)
set_property (GLOBAL PROPERTY USE_FOLDERS ON)


find_package(VTK REQUIRED)
find_package(ITK REQUIRED)


set(CMAKE_MODULE_PATH
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake
  ${Slicer_CMAKE_DIR}
  ${CMAKE_MODULE_PATH}
  )

# build options
if(WIN32)
  option(BUILD_CONSOLE "Build ${PROJECT_NAME} executable as a console app on windows (allows debug output)" ON)
endif(WIN32)

# set libs
set(VISUAL_LIBS 
  vtkDMRI
  vtkTeem
)


# set apps
set(VISUAL_APPS )

if(BUILD_VisualNeurosurgery)
  set(VISUAL_APPS VisualNeurosurgery ${VISUAL_APPS} )
endif()


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# add libs directory
foreach(lib ${VISUAL_LIBS})
  set(libExportStr "VISUAL_${lib}_EXPORT")
  string(TOUPPER ${libExportStr} MY_LIBRARY_EXPORT_DIRECTIVE)
  set(MY_EXPORT_HEADER_PREFIX ${lib})
  string(TOUPPER ${lib} MY_LIBNAME)
  configure_file(
    ${CMAKE_SOURCE_DIR}/Libs/visualSDKExport.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/Libs/${lib}/${lib}Export.h
    )
  add_subdirectory(Libs/${lib})
endforeach()

# add apps directory
foreach(app ${VISUAL_APPS})
  add_subdirectory(Applications/${app})
endforeach()

#-----------------------------------------------------------------------------
# The commands in this directory are intended to be executed as
# the end of the whole configuration process, as a "last step".
# This directory is typically the last add_subdirectory in the main CMakeLists.txt.
#add_subdirectory(CMake/LastConfigureStep)
