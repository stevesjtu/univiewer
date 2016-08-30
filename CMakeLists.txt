cmake_minimum_required(VERSION 2.8)
 
PROJECT(univiewer)
set(VTK_DIR /Users/rebecca/SteveShi/Library/VTK/lib/cmake/vtk-7.0)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(/Users/rebecca/SteveShi/Library)
 
add_executable(univiewer univiewer)
 
target_link_libraries(univiewer ${VTK_LIBRARIES})
