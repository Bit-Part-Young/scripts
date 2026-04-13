#!/bin/bash

# LAMMPS 2026 版本的会报错

module purge

# SiYuan
# module load intel-oneapi-compilers/2021.4.0
# module load intel-oneapi-mkl/2021.4.0
# module load intel-oneapi-mpi/2021.12.1
# module load intel-oneapi-tbb/2021.4.0

cmake -C ../cmake/presets/all_on.cmake -C ../cmake/presets/nolib.cmake \
      -D CMAKE_C_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_COMPILER_LAUNCHER=ccache \
      -D FFT=MKL \
      -D BUILD_SHARED_LIBS=yes \
      -C ../cmake/presets/intel.cmake ../cmake

cmake -D PKG_ML-PACE=on -D LOCAL_ML-PACE="${HOME}/opt/lammps-user-pace" ../cmake
# cmake -D PKG_ML-PACE=on -D LOCAL_ML-PACE="${HOME}/yangsl/opt/lammps-user-pace" ../cmake

# cmake -D PKG_ML-QUIP=yes ../cmake
cmake -D PKG_ML-QUIP=yes -D DOWNLOAD_QUIP=no -D QUIP_LIBRARY=${HOME}/opt/quip-build-ifort-icc/src/quip_build/build/lammps/libquip.a ../cmake

cmake -D PKG_USER-NEP=on ../cmake
