project(DWIConvert)

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
# set(configure_header_file vtkDMRIConfigure.h)
# configure_file(
#   ${CMAKE_CURRENT_SOURCE_DIR}/${configure_header_file}.in
#   ${CMAKE_CURRENT_BINARY_DIR}/${configure_header_file}
#   )

# --------------------------------------------------------------------------
# Dependencies
# --------------------------------------------------------------------------

#
#ITK
#
#include(${ITK_USE_FILE})


# --------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------
set(${PROJECT_NAME}_SRCS
  #FSLToNrrd.cxx
  #NrrdToFSL.cxx
  NRRDDWIConverter.h
  NRRDDWIConverter.cxx
  FSLDWIConverter.h
  FSLDWIConverter.cxx
  DWIConvertUtils.h
  DWIConvertUtils.cxx
  DWIConverterFactory.h
  DWIConverterFactory.cxx
  DWIConverter.h
  DWIConverter.cxx
  DWIDICOMConverterBase.h
  DWIDICOMConverterBase.cxx
  GEDWIConverter.h
  GenericDWIConverter.h
  HitachiDWIConverter.h
  PhilipsDWIConverter.h
  SiemensDWIConverter.h
  StringContains.h
  GEDWIConverter.cxx
  GenericDWIConverter.cxx
  HitachiDWIConverter.cxx
  PhilipsDWIConverter.cxx
  SiemensDWIConverter.cxx
  SiemensDWIConverter.hxx
  DWIConvertLib.cxx
  DWIConvertLib.h
  DWIMetaDataDictionaryValidator.h
  DWIMetaDataDictionaryValidator.cxx
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

set(srcs ${${PROJECT_NAME}_SRCS})

add_library(${lib_name} STATIC  ${srcs})



target_include_directories(${PROJECT_NAME}
    PUBLIC
    $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    PRIVATE
    
    )

set(libs
  ${ITK_LIBRARIES}
  )
target_link_libraries(${lib_name} ${libs})


# --------------------------------------------------------------------------
# Folder
# --------------------------------------------------------------------------
if(NOT DEFINED ${lib_name}_FOLDER)
  set(${lib_name}_FOLDER "Libs/DWIConver")
endif()
if(NOT "${${lib_name}_FOLDER}" STREQUAL "")
  set_target_properties(${lib_name} PROPERTIES FOLDER ${${lib_name}_FOLDER})
endif()


# --------------------------------------------------------------------------
# exe
# --------------------------------------------------------------------------
set(TEST_FOLDER
Libs/DWIConver/Testing
)
set(PROGS_LIST_1
  Convert4DImageTo3DSeries
  ConvertDicomToDWI
  )
foreach(prog ${PROGS_LIST_1})
  add_executable(${prog} ${prog}.cxx)
  target_link_libraries(${prog} ${lib_name})
  set_target_properties(${prog} PROPERTIES FOLDER ${TEST_FOLDER})
endforeach()


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

