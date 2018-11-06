cmake_minimum_required(VERSION 2.8)
 
PROJECT(univiewer)
IF (WIN32)
  MESSAGE(STATUS "Configure on Microsoft Windows.")
  set(VTK_DIR "C:/Users/CNJISHI10/Development/VTK/lib/cmake/vtk-8.1")
ELSEIF (APPLE)
  MESSAGE(STATUS "Configure on Apple macOS.")
  set(VTK_DIR "/Users/rebecca/SteveShi/Library/VTK/lib/cmake/vtk-7.0")
ELSEIF (UNIX)
	MESSAGE(STATUS "Configure on UNIX-like OS.")
ENDIF ()

set(CMAKE_CXX_STANDARD 14)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories("C:/Users/CNJISHI10/Development/VTK/include/vtk-8.1")

aux_source_directory(./ ALL_SOURCE)

add_executable(univiewer ${ALL_SOURCE})
 
target_link_libraries(univiewer ${VTK_LIBRARIES})
