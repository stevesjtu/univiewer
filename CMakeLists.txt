cmake_minimum_required(VERSION 2.8)
 
PROJECT(univiewer)
set(VTK_DIR "C:/Users/CNJISHI10/Development/VTK/lib/cmake/vtk-8.1")

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories("C:/Users/CNJISHI10/Development/VTK/include/vtk-8.1")

aux_source_directory(./ ALL_SOURCE)

add_executable(univiewer ${ALL_SOURCE})
 
target_link_libraries(univiewer ${VTK_LIBRARIES})
