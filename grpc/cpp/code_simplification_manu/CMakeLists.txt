# Template Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( code_simplification_manu )

cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )
find_package(Boost REQUIRED COMPONENTS system program_options)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

project(TEST C CXX)
add_executable(TEST test.cpp)
target_link_libraries(TEST Boost::system Boost::program_options)


