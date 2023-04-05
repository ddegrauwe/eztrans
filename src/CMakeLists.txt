file( GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/*.F90 )

# define library
add_library( eztrans STATIC ${sources} )

# install library
install( TARGETS ${PROJECT_NAME} DESTINATION lib/ )

# install modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include )
