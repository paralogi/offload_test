# OFFLOAD_TEST

Проект для тестирования способов ускорения вычислений и новых стандартов программирования.

# Структура проекта

* **src/openmp**

Тест `OpenMP` для `NVPTX`-`Cuda`.

* **src/sycl**

Тест `Sycl` для `OpenCL` / `Cuda`.

* **src/eigen3**

Тест `Eigen3` совместно с `Intel Math Kernel Library (oneMKL)`.

# Настройка среды

Системные требования:
* [gcc 7.1.0](https://gcc.gnu.org/install/)
* [cmake 3.14](https://gcc.gnu.org/install/)
* [cuda 11.6](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)

Установка зависимостей:

```sh
sudo apt install -y build-essential cmake libelf-dev libffi-dev pkg-config ninja-build git
```

Сборка компилятора:

```sh
git clone -b sycl https://github.com/intel/llvm llvm-intel-sycl
mkdir llvm-intel-sycl/build
cd llvm-intel-sycl/build
```

```sh
export CUDA_LIB_PATH=/usr/local/cuda-11.6/lib64/stubs
export ROOT=$(pwd)/..
cmake ${ROOT}/llvm -GNinja \
    -DCUDA_TOOLKIT_ROOT_DIR="/usr/local/cuda-11.6" \
    -DCMAKE_INSTALL_PREFIX="/opt/llvm-intel-sycl" \
    -DLLVM_EXTERNAL_SYCL_SOURCE_DIR="${ROOT}/sycl" \
    -DLLVM_EXTERNAL_LLVM_SPIRV_SOURCE_DIR="${ROOT}/llvm-spirv" \
    -DLLVM_EXTERNAL_LIBDEVICE_SOURCE_DIR="${ROOT}/libdevice" \
    -DCMAKE_BUILD_TYPE="RELEASE" \
    -DLLVM_ENABLE_PROJECTS="clang;openmp;sycl;llvm-spirv;opencl;libdevice;libclc" \
    -DLLVM_EXTERNAL_PROJECTS="sycl;llvm-spirv;opencl;libdevice" \
    -DLLVM_TARGETS_TO_BUILD="X86;NVPTX" \
    -DLIBCLC_TARGETS_TO_BUILD="nvptx64--;nvptx64--nvidiacl" \
    -DCLANG_OPENMP_NVPTX_DEFAULT_ARCH="sm_60" \
    -DLIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES="37,60,70" \
    -DLIBCLC_GENERATE_REMANGLED_VARIANTS=ON \
    -DSYCL_BUILD_PI_CUDA=ON \
    -DLLVM_BUILD_TOOLS=ON \
    -DLLVM_ENABLE_ASSERTIONS=OFF \
    -DLLVM_ENABLE_DUMP=OFF \
    -DLLVM_INCLUDE_TESTS=OFF \
    -DLLVM_INCLUDE_DOCS=OFF 
```
    
```sh
make -j$(nproc)
sudo make -j$(nproc) install 
```

Переменные среды:

```sh
CUDA_PATH=/usr/local/cuda-11.6
LLVM_PATH=/opt/llvm-intel-sycl

CMAKE_PREFIX_PATH="$CUDA_PATH/lib64:$CMAKE_PREFIX_PATH"
CMAKE_PREFIX_PATH="$LLVM_PATH/lib:$CMAKE_PREFIX_PATH"
export CMAKE_PREFIX_PATH

LD_LIBRARY_PATH="$CUDA_PATH/lib64:$CUDA_PATH/lib64/stubs:$LD_LIBRARY_PATH"
LD_LIBRARY_PATH="$LLVM_PATH/lib:$LLVM_PATH/libexec:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH

PATH="$CUDA_PATH/bin:$PATH"
PATH="$LLVM_PATH/bin:$PATH"
export PATH
```

# Ссылки

* [OpenMP 4.5 Examples](https://openmp.org/wp-content/uploads/openmp-examples-4.5.0.pdf)
* [Offloading Support in GCC](https://gcc.gnu.org/wiki/Offloading)
* [LLVM/OpenMP Documentation and Support](https://openmp.llvm.org/SupportAndFAQ.html)
* [Building LLVM/Clang with OpenMP Offloading to NVIDIA GPUs](https://hpc-wiki.info/hpc/Building_LLVM/Clang_with_OpenMP_Offloading_to_NVIDIA_GPUs)
* [How To Build And Run Your Modern Parallel Code In C++17 and OpenMP 4.5 Library On NVIDIA GPUs](https://devmesh.intel.com/posts/724749/how-to-build-and-run-your-modern-parallel-code-in-c-17-and-openmp-4-5-library-on-nvidia-gpus)
