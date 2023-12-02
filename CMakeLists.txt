cmake_minimum_required( VERSION 3.5 )
project( OFFLOAD_TEST LANGUAGES CXX C )

add_subdirectory( src/openmp )
add_subdirectory( src/sycl )
#add_subdirectory( src/eigen3 )
