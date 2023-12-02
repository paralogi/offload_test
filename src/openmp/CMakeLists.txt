set( CLANG_TOOLCHAIN True )

if ( CLANG_TOOLCHAIN )
    set( CMAKE_CXX_COMPILER clang++ )
    set( OPENMP_OPTIONS -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda )
else()
    set( CMAKE_CXX_COMPILER g++ )
    set( OPENMP_OPTIONS -fopenmp -foffload=nvptx-none -foffload=-misa=sm_35 -fno-stack-protector -fcf-protection=none )
endif()

add_executable( offload_openmp main.cpp )
target_compile_options( offload_openmp PRIVATE ${OPENMP_OPTIONS} )
target_link_options( offload_openmp PRIVATE ${OPENMP_OPTIONS} )
