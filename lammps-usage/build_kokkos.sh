#!/bin/bash

# LAMMPS 2026 版本的 Master 上无法运行（调用 1、2 块 GPU 均报错），在 node2 节点正常运行

module purge
module load OpenMPI/5.0.10

# 不构建动态库
cmake -C ../cmake/presets/all_on.cmake -C ../cmake/presets/nolib.cmake \
      -D PKG_INTEL=no \
      -C ../cmake/presets/kokkos-cuda-nowrapper.cmake ../cmake

cmake -D PKG_ML-PACE=on -D LOCAL_ML-PACE="${HOME}/opt/lammps-user-pace" ../cmake
