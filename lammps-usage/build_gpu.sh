#!/bin/bash

module purge
module load OpenMPI/5.0.10

# 不构建动态库
cmake -C ../cmake/presets/all_on.cmake -C ../cmake/presets/nolib.cmake \
      -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache \
      -D PKG_INTEL=no \
      -C ../cmake/presets/gpu-cuda.cmake ../cmake
