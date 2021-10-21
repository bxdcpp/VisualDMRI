project(VisualDMRITesting)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/Libs
  ${CMAKE_SOURCE_DIR}/Libs/vtkDMRI
)

set(header_files
)

set(source_files
  main.cpp
)


add_executable(${PROJECT_NAME} WIN32
  ${header_files}
  ${source_files}
)

# if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
#   if(NOT show_console_debug)
#     set_property(TARGET ${PROJECT_NAME} PROPERTY WIN32_EXECUTABLE true)
#   endif()
# endif()

#Set_Target_Properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#Set_Target_Properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS_DEBUG "/SUBSYSTEM:CONSOLE /ENTRY:mainCRTStartup")

target_link_libraries(${PROJECT_NAME}
    vtkDMRI
)

