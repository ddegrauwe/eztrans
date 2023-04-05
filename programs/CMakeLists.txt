# define executable
add_executable( driver driver.F90 )

message("CMAKE_Fortran_MODULE_DIRECTORY = ${CMAKE_Fortran_MODULE_DIRECTORY}")

target_include_directories(driver PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY}) 
target_link_libraries(driver eztrans)

# install library
install( TARGETS driver DESTINATION bin/ )
