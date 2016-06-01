cmake_minimum_required(VERSION 2.8)
 
PROJECT(univiewer)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(univiewer univiewer)
 
target_link_libraries(univiewer ${VTK_LIBRARIES})