find_package( MKL REQUIRED CONFIG )
find_package( Eigen3 REQUIRED CONFIG )

add_executable( offload_eigen3 main.cpp )
target_include_directories( offload_eigen3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries( offload_eigen3 PRIVATE MKL::MKL Eigen3::Eigen )
target_compile_definitions( offload_eigen3 PRIVATE -DEIGEN_USE_MKL_ALL )
