set( CMAKE_CXX_COMPILER clang++ )
set( SYCL_OPTIONS -fsycl -fsycl-targets=nvptx64-nvidia-cuda -fsycl-unnamed-lambda -Wno-linker-warnings )

add_executable( offload_sycl main.cpp )
target_compile_options( offload_sycl PRIVATE ${SYCL_OPTIONS} )
target_link_options( offload_sycl PRIVATE ${SYCL_OPTIONS} )
