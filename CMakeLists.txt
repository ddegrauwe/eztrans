cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

project( eztrans LANGUAGES Fortran)

# module directory
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include/ )

# sources
add_subdirectory(src)

# programs
add_subdirectory(programs)
