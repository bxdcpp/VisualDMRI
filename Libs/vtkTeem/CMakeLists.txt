project(vtkTeem)

#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.13.4)
#-----------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------
if(NOT DEFINED BUILD_SHARED_LIBS)
  option(BUILD_SHARED_LIBS "Build with shared libraries." ON)
endif()

# --------------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------------

#
# VTK
#
#find_package(VTK REQUIRED)
if(${VTK_VERSION} VERSION_LESS "8.90.0")
  include(${VTK_USE_FILE})
endif()

#
# ITK
#

#set(${PROJECT_NAME}_ITK_COMPONENTS
#  ITKCommon
#  ITKVNL
#  )
#find_package(ITK 5.0 COMPONENTS ${${PROJECT_NAME}_ITK_COMPONENTS} REQUIRED)
include(${ITK_USE_FILE})

#
# Teem
#
find_package(Teem REQUIRED NO_MODULE)
# XXX Since no TeemConfigVersion.cmake is configured, let's check the version explicitly
if(${Teem_VERSION_MAJOR}.${Teem_VERSION_MINOR}.${Teem_VERSION_PATCH} VERSION_LESS 1.10)
  message(FATAL_ERROR
    "Could not find a configuration file for package \"Teem\" that is compatible with requested version \"1.10\".")
endif()
include(${Teem_USE_FILE})

# XXX Teem 1.11.0 or Teem >= 4168 - This variable won't be needed anymore.
#     See http://sourceforge.net/p/teem/code/4168/
set(Teem_LIBRARIES teem)

# --------------------------------------------------------------------------
# Configure headers
# --------------------------------------------------------------------------
set(configure_header_file vtkTeemConfigure.h)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/${configure_header_file}.in
  ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
  )

# --------------------------------------------------------------------------
# Install headers
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  set(${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT ON)
endif()
if(NOT ${PROJECT_NAME}_INSTALL_NO_DEVELOPMENT)
  file(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/*.(h|txx)")
  install(
    FILES ${headers} ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
    DESTINATION include/${PROJECT_NAME} COMPONENT Development)
endif()

# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------
set(vtkTeem_SRCS
  vtkDiffusionTensorMathematics.cxx
  vtkDiffusionTensorGlyph.cxx
  vtkTeemNRRDReader.cxx
  vtkTeemNRRDWriter.cxx
  vtkImageLabelCombine.cxx
  )

# --------------------------------------------------------------------------
# Include dirs
# --------------------------------------------------------------------------
#set(include_dirs
#  ${CMAKE_CURRENT_SOURCE_DIR}
#  ${CMAKE_CURRENT_BINARY_DIR}
#  )
#include_directories(BEFORE ${include_dirs})


# --------------------------------------------------------------------------
# Build library
# --------------------------------------------------------------------------
set(lib_name ${PROJECT_NAME})

set(srcs ${vtkTeem_SRCS})
add_library(${lib_name} ${srcs})

set(libs
  #itkvnl
  ITKCommon
  ${Teem_LIBRARIES}
  ${VTK_LIBRARIES}
  )
target_link_libraries(${lib_name} ${libs})

target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${Teem_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${VTK_INCLUDE_DIRS}> 
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    
    )



# Apply user-defined properties to the library target.
#if(Slicer_LIBRARY_PROPERTIES)
#  set_target_properties(${lib_name} PROPERTIES ${Slicer_LIBRARY_PROPERTIES})
#endif()

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
if(NOT DEFINED ${PROJECT_NAME}_EXPORT_FILE)
  set(${PROJECT_NAME}_EXPORT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)
endif()
#export(TARGETS ${lib_name} APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})
export(TARGETS ${lib_name} FILE ${${PROJECT_NAME}_EXPORT_FILE})

# --------------------------------------------------------------------------
# Install library
# --------------------------------------------------------------------------
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_BIN_DIR)
  set(${PROJECT_NAME}_INSTALL_BIN_DIR bin)
endif()
if(NOT DEFINED ${PROJECT_NAME}_INSTALL_LIB_DIR)
  set(${PROJECT_NAME}_INSTALL_LIB_DIR lib/${PROJECT_NAME})
endif()

install(TARGETS ${lib_name}
  RUNTIME DESTINATION ${${PROJECT_NAME}_INSTALL_BIN_DIR} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${${PROJECT_NAME}_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${${PROJECT_NAME}_INSTALL_LIB_DIR} COMPONENT Development
  )

# --------------------------------------------------------------------------
# Python Wrapping
# --------------------------------------------------------------------------
# if(VTK_WRAP_PYTHON)
#   include(vtkMacroKitPythonWrap)

#   vtkMacroKitPythonWrap(
#     KIT_NAME ${lib_name}
#     KIT_SRCS ${vtkTeem_SRCS}
#     KIT_INSTALL_BIN_DIR ${${PROJECT_NAME}_INSTALL_BIN_DIR}
#     KIT_INSTALL_LIB_DIR ${${PROJECT_NAME}_INSTALL_LIB_DIR}
#     )
#   # Export target
#   export(TARGETS ${lib_name}Python ${lib_name}PythonD APPEND FILE ${${PROJECT_NAME}_EXPORT_FILE})
#   # Folder
#   if(NOT "${${PROJECT_NAME}_FOLDER}" STREQUAL "")
#     set_target_properties(${lib_name}Python PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
#     set_target_properties(${lib_name}PythonD PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
#     if(TARGET ${lib_name}Hierarchy)
#       set_target_properties(${lib_name}Hierarchy PROPERTIES FOLDER ${${PROJECT_NAME}_FOLDER})
#     endif()
#   endif()
# endif()

# --------------------------------------------------------------------------
# Testing (requires some of the examples)
# --------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(Testing)
endif()

# --------------------------------------------------------------------------
# Set INCLUDE_DIRS variable
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
  CACHE INTERNAL "${PROJECT_NAME} include dirs" FORCE)

# --------------------------------------------------------------------------
# Configure support file
# --------------------------------------------------------------------------
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/vtkTeemConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/vtkTeemConfig.cmake
  )
