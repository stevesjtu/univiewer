cmake_minimum_required(VERSION 3.12)
 
PROJECT(univiewer)
# IF (WIN32)
#   MESSAGE(STATUS "Configure on Microsoft Windows.")
#   set(VTK_DIR "C:/Users/CNJISHI10/Development/VTK/lib/cmake/vtk-8.1")
# ELSEIF (APPLE)
#   MESSAGE(STATUS "Configure on Apple macOS.")
#   set(VTK_DIR "/Users/rebecca/SteveShi/Development/vtk-8.1/lib/cmake/vtk-8.1")
# ELSEIF (UNIX)
# 	MESSAGE(STATUS "Configure on UNIX-like OS.")
# ENDIF ()

set(VTK_DIR "D:/repos/common/vtk/lib/cmake/vtk-9.0")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -m64 -mavx2")

find_package(VTK REQUIRED)

set(HDF5_USE_STATIC_LIBRARIES OFF)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
add_definitions(-D H5_BUILT_AS_DYNAMIC_LIB)

aux_source_directory(./ ALL_SOURCE)

MESSAGE(STATUS "Build type is ${CMAKE_BUILD_TYPE}.")
MESSAGE(STATUS "##########################################")

if( ${CMAKE_BUILD_TYPE} MATCHES Release )
  add_executable(univiewer ${ALL_SOURCE})
elseif (${CMAKE_BUILD_TYPE} MATCHES Debug)
  add_executable(univiewer ${ALL_SOURCE})
endif()

add_executable(univiewer ${ALL_SOURCE})
target_link_libraries(univiewer ${VTK_LIBRARIES} ${HDF5_CXX_LIBRARIES})
vtk_module_autoinit(
  TARGETS univiewer
  MODULES ${VTK_LIBRARIES}
  )